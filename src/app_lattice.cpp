/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "spktype.h"
#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "domain.h"
#include "lattice.h"
#include "random_mars.h"
#include "random_park.h"
#include "cluster.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define DELTA 10000

enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};

/* ---------------------------------------------------------------------- */

AppLattice::AppLattice(SPPARKS *spk, int narg, char **arg) : App(spk,narg,arg)
{
  appclass = LATTICE;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  sectorflag = 0;
  nset = 0;
  set = NULL;
  nstop = 1.0;
  tstop = 0.0;

  sweepflag = NOSWEEP;
  ranapp = NULL;
  ranstrict = NULL;
  siteseeds = NULL;
  sitelist = NULL;
  Lmask = false;
  mask = NULL;

  allow_update = 0;

  temperature = 0.0;

  propensity = NULL;
  i2site = NULL;

  comm = NULL;
  sweep = NULL;

  nlocal = nghost = nmax = 0;
  owner = NULL;
  index = NULL;

  maxneigh = 0;
  numneigh = NULL;
  neighbor = NULL;

  dt_sweep = 0.0;
  naccept = nattempt = 0;
  nsweeps = 0;
  
  update_only = 0;
}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  Solve *s;
  for (int i = 0; i < nset; i++) {
    s = free_set(i);
    delete s;
  }
  delete [] set;

  delete ranapp;
  delete ranstrict;
  memory->destroy(siteseeds);
  memory->destroy(sitelist);
  memory->destroy(mask);

  delete comm;

  memory->destroy(owner);
  memory->destroy(index);

  memory->destroy(numneigh);
  memory->destroy(neighbor);
}

/* ---------------------------------------------------------------------- */

void AppLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"sweep") == 0) set_sweep(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"update_only") == 0) set_update_only(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  // error checks

  if (solve == NULL && sweepflag == NOSWEEP)
    error->all(FLERR,"App needs a KMC or rejection KMC solver");
  if (solve && sweepflag != NOSWEEP)
    error->all(FLERR,"App cannot use both a KMC and rejection KMC solver");

  if (solve && allow_kmc == 0)
    error->all(FLERR,"KMC events are not implemented in app");
  if (sweepflag != NOSWEEP && allow_rejection == 0)
    error->all(FLERR,"Rejection events are not implemented in app");
  if (sweepflag != NOSWEEP && Lmask && allow_masking == 0)
    error->all(FLERR,"Mask logic not implemented in app");

  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all(FLERR,"Cannot use KMC solver in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
    error->all(FLERR,"Cannot use random rejection KMC in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
    error->all(FLERR,"Cannot use raster rejection KMC in parallel with no sectors");
  if (sectorflag && sweepflag == COLOR_STRICT)
    error->all(FLERR,"Cannot use color/strict rejection KMC with sectors");

  if (sweepflag && dt_sweep == 0.0) error->all(FLERR,"App did not set dt_sweep");

  // if sectors, set number of sectors

  int dimension = domain->dimension;

  nsector = 1;
  if (sectorflag) {
    if (nsector_user) nsector = nsector_user;
    else if (dimension == 1) nsector = 2;
    else if (dimension == 2) nsector = 4;
    else nsector = 8;

    if (dimension == 3) {
      if (nsector == 2 && (domain->procgrid[1] != 1 || 
			   domain->procgrid[2] != 1))
	error->all(FLERR,"Invalid number of sectors");
      if (nsector == 4 && domain->procgrid[2] != 1)
	error->all(FLERR,"Invalid number of sectors");
    }
    if (dimension == 2) {
      if (nsector == 2 && domain->procgrid[1] != 1)
	error->all(FLERR,"Invalid number of sectors");
      if (nsector == 8)
	error->all(FLERR,"Invalid number of sectors");
    }
    if (dimension == 1 && nsector != 2)
      error->all(FLERR,"Invalid number of sectors");
  }

  // if coloring, determine number of colors
  // setup test for create_set
  // check periodicity against lattice extent

  ncolors = 1;
  if (sweepflag == COLOR || sweepflag == COLOR_STRICT) {
    int delcolor = delevent + delpropensity;
    if (domain->lattice == NULL)
      error->all(FLERR,"Cannot color without a lattice definition of sites");
    if (contiguous_sites() == 0)
      error->all(FLERR,"Cannot color without contiguous site IDs");
    ncolors = domain->lattice->ncolors(delcolor);
    if (ncolors == 0)
      error->all(FLERR,"Cannot color this combination of lattice and app");
  }

  if (nsector > 1 && ncolors > 1) bothflag = 1;
  else bothflag = 0;

  // create sets based on sectors and coloring
  // set are either all sectors or all colors or both
  // for both, first sets are entire sections, remaining are colors in sectors
  // if new nset is same as old nset, pass each set's solver to create_set,
  //   so it can reuse solver and its RNG,
  //   so consecutive runs are identical to one continuous run

  int nsetold = nset;
  Solve **sold = new Solve*[nsetold];

  for (int i = 0; i < nset; i++) sold[i] = free_set(i);
  delete [] set;

  if (nsector == 1 && ncolors == 1) {
    nset = 1;
    set = new Set[nset];
    if (nset == nsetold) create_set(0,0,0,sold[0]);
    else create_set(0,0,0,NULL);
  } else if (nsector > 1 && ncolors == 1) {
    nset = nsector;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) {
      if (nset == nsetold) create_set(i,i+1,0,sold[i]);
      else create_set(i,i+1,0,NULL);
    }
  } else if (ncolors > 1 && nsector == 1) {
    nset = ncolors;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) {
      if (nset == nsetold) create_set(i,0,i+1,sold[i]);
      else create_set(i,0,i+1,NULL);
    }
  } else if (bothflag) {
    nset = nsector + ncolors*nsector;
    set = new Set[nset];
    int m = 0;
    for (int i = 0; i < nsector; i++) {
      if (nset == nsetold) create_set(m,i+1,0,sold[m]);
      else create_set(m,i+1,0,NULL);
      m++;
    }
    for (int i = 0; i < nsector; i++) 
      for (int j = 0; j < ncolors; j++) {
	if (nset == nsetold) create_set(m,i+1,j+1,sold[m]);
        else create_set(m,i+1,j+1,NULL);
	m++;
      }
  }

  if (nset != nsetold)
    for (int i = 0; i < nsetold; i++) delete sold[i];
  delete [] sold;

  // initialize mask array

  if (!Lmask && mask) {
    memory->destroy(mask);
    mask = NULL;
  }
  if (Lmask && mask == NULL) {
    memory->create(mask,nlocal+nghost,"app:mask");
    for (int i = 0; i < nlocal+nghost; i++) mask[i] = 0;
  }

  // setup RN generators, only on first init
  // ranapp is used for all options except sweep color/strict
  // setup ranapp so different on every proc
  // if color/strict, initialize per-lattice site seeds

  if (ranapp == NULL) {
    ranapp = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    ranapp->reset(seed,me,100);
  }

  if (sweepflag != COLOR_STRICT) {
    delete ranstrict;
    memory->destroy(siteseeds);
    ranstrict = NULL;
    siteseeds = NULL;
  }

  if (sweepflag == COLOR_STRICT && ranstrict == NULL) {
    ranstrict = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    memory->create(siteseeds,nlocal,"app:siteseeds");
    for (int i = 0; i < nlocal; i++) {
      ranstrict->tagreset(seed,id[i],100);
      siteseeds[i] = ranstrict->seed;
    }
  }

  // initialize comm, both for this proc's full domain and sectors
  // recall comm->init in case sectoring has changed

  if (comm == NULL) comm = new CommLattice(spk);
  comm->init(nsector,delpropensity,delevent,NULL);

  // set sweep function ptr

  if (sweepflag != NOSWEEP) {
    if (sweepflag != COLOR_STRICT && !Lmask)
      sweep = &AppLattice::sweep_nomask_nostrict;
    else if (sweepflag != COLOR_STRICT && Lmask)
      sweep = &AppLattice::sweep_mask_nostrict;
    else if (sweepflag == COLOR_STRICT && !Lmask)
      sweep = &AppLattice::sweep_nomask_strict;
    else if (sweepflag == COLOR_STRICT && Lmask)
      sweep = &AppLattice::sweep_mask_strict;
  } else sweep = NULL;

  // app-specific initialization, after general initialization

  init_app();

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice::setup()
{
  // app-specific setup, before propensities are computed

  setup_app();

  // initialize propensities for KMC solver within each set
  // comm insures ghost sites are up to date

  if (solve) {
    comm->all();
    for (int i = 0; i < nset; i++) {
      for (int m = 0; m < set[i].nlocal; m++)
	set[i].propensity[m] = site_propensity(set[i].site2i[m]);
      set[i].solve->init(set[i].nlocal,set[i].propensity);
    }
  }

  // convert per-sector time increment info to KMC params

  if (sectorflag && solve) {
    if (tstop > 0.0) {
      Ladapt = false;
      dt_kmc = tstop;
    }

    if (nstop > 0.0) {
      Ladapt = true;
      double pmax = 0.0;
      for (int i = 0; i < nset; i++) {
	int ntmp = set[i].solve->get_num_active();
	if (ntmp > 0) {
	  double ptmp = set[i].solve->get_total_propensity();
	  ptmp /= ntmp;
	  pmax = MAX(ptmp,pmax);
	}
      }
      double pmaxall;
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
    }

    dt_kmc = MIN(dt_kmc,stoptime-time);
  }

  // convert rejection info to rKMC params
  // nloop and nselect are set whether sectoring is used or not
  // if bothflag, list of active sets starts with nsector

  if (sweepflag != NOSWEEP) {
    int first = 0;
    if (bothflag) first = nsector;

    if (sweepflag == RANDOM) {
      if (nstop > 0.0) {
	for (int i = first; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (nstop*set[i].nlocal);
	}
      }
      if (tstop > 0.0) {
	double n = tstop / (dt_sweep/nglobal);
	for (int i = first; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (n/nglobal * set[i].nlocal);
	}
      }

    } else if (sweepflag == RASTER || 
	       sweepflag == COLOR || sweepflag == COLOR_STRICT) {
      int n;
      if (nstop > 0.0) n = static_cast<int> (nstop);
      if (tstop > 0.0) n = static_cast<int> (tstop/dt_sweep);
      for (int i = first; i < nset; i++) {
	set[i].nloop = n;
	set[i].nselect = n * set[i].nlocal;
      }
    }

    double nme = 0.0;
    for (int i = first; i < nset; i++) nme += set[i].nselect;
    double ntotal;
    MPI_Allreduce(&nme,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);

    dt_rkmc = ntotal/nglobal * dt_sweep;
    if (dt_rkmc == 0.0)
      error->all(FLERR,"Choice of sector stop led to no rKMC events");
    dt_rkmc = MIN(dt_rkmc,stoptime-time);
  }

  // setup sitelist if sweepflag = RANDOM
  // do this every run since sector timestep could have changed

  if (sweepflag == RANDOM) {
    memory->destroy(sitelist);
    int n = 0;
    for (int i = 0; i < nset; i++) n = MAX(n,set[i].nselect);
    memory->create(sitelist,n,"app:sitelist");
  }

  // second stage of app-specific setup

  setup_end_app();

  // setup future output

  nextoutput = output->setup(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice::iterate()
{
  timer->barrier_start(TIME_LOOP);

  if (solve) {
    if (sectorflag == 0) 
      iterate_kmc_global(stoptime);
    else if (allow_update && update_only)
      iterate_update_only(stoptime,dt_kmc);
    else
      iterate_kmc_sector(stoptime);
  } else {
    if (allow_update && update_only)
      iterate_update_only(stoptime,dt_rkmc);
    else
      iterate_rejection(stoptime);
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   KMC solver on entire domain
   can only be invoked in serial
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_global(double stoptime)
{
  int isite;
  
  // global KMC runs with one set
  // save ptr to system solver

  Solve *hold_solve = solve;
  solve = set[0].solve;
  propensity = set[0].propensity;
  i2site = set[0].i2site;

  int done = 0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt_step);
    timer->stamp(TIME_SOLVE);

    if (isite >= 0) {
      time += dt_step;
      if (time <= stoptime) {
	site_event(isite,ranapp);
	naccept++;
	timer->stamp(TIME_APP);
      } else {
	done = 1;
	time = stoptime;
      }
    } else {
      done = 1;
      time = stoptime;
    }

    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }

  // restore system solver

  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   KMC solver on sectors
   can be invoked in serial or parallel
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_sector(double stoptime)
{
  int i,isite,done;
  double dt,timesector;
  double pmax,pmaxall;

  // save ptr to system solver

  Solve *hold_solve = solve;

  int alldone = 0;
  while (!alldone) {
    if (Ladapt) pmax = 0.0;

    for (int iset = 0; iset < nset; iset++) {
      timer->stamp();

      if (nprocs > 1) {
	comm->sector(iset);
	timer->stamp(TIME_COMM);
      }

      solve = set[iset].solve;
      
      propensity = set[iset].propensity;
      i2site = set[iset].i2site;
      int *site2i = set[iset].site2i;

      // update propensities for sites which neighbor a site outside sector
      // necessary since outside sites may have changed
      // attribute this chunk of time to comm, b/c due to decomposition

      int *bsites = set[iset].bsites;
      int *border = set[iset].border;
      int nborder = set[iset].nborder;

      int nsites = 0;
      for (int m = 0; m < nborder; m++) {
	i = border[m];
	isite = i2site[i];
	bsites[nsites++] = isite;
	propensity[isite] = site_propensity(i);
      }
      
      solve->update(nsites,bsites,propensity);
      timer->stamp(TIME_COMM);
      
      // pmax = maximum sector propensity per site

      if (Ladapt) {
	int ntmp = solve->get_num_active();
	if (ntmp > 0) {
	  double ptmp = solve->get_total_propensity();
	  ptmp /= ntmp;
	  pmax = MAX(ptmp,pmax);
	}
      }
      
      // execute events until sector time threshhold reached
      
      done = 0;
      timesector = 0.0;
      while (!done) {
	timer->stamp();
	isite = solve->event(&dt);
	timer->stamp(TIME_SOLVE);
	
	if (isite < 0) done = 1;
	else {
	  timesector += dt;
	  if (timesector >= dt_kmc) done = 1;
	  else {
	    site_event(site2i[isite],ranapp);
	    naccept++;
	  }
	  timer->stamp(TIME_APP);
	}
      }
      
      if (nprocs > 1) {
	comm->reverse_sector(iset);
	timer->stamp(TIME_COMM);
      }
    }

    if (allow_update) user_update(dt_kmc);

    // keep looping until overall time threshhold reached

    nsweeps++;
    time += dt_kmc;
    if (time >= stoptime) alldone = 1;
    if (alldone || time >= nextoutput)
      nextoutput = output->compute(time,alldone);
    timer->stamp(TIME_OUTPUT);

    // recompute dt_kmc if adaptive, based on pmax across all sectors

    if (Ladapt) {
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
      dt_kmc = MIN(dt_kmc,stoptime-time);
    }
  }

  // restore system solver

  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   rejection KMC solver
 ------------------------------------------------------------------------- */

void AppLattice::iterate_rejection(double stoptime)
{
  int i,icolor,nselect,nrange,jset;
  int *site2i;

  // set loop is over:
  // sectors if there are sectors and no colors
  // colors if there are colors and no sectors
  // first nsector sets if there are both sectors and colors

  int nset_loop = nset;
  if (bothflag) nset_loop = nsector;

  int done = 0;
  while (!done) {
    for (int iset = 0; iset < nset_loop; iset++) {
      if (nprocs > 1) {
	timer->stamp();
	if (sectorflag) comm->sector(iset);
	else comm->all();
	timer->stamp(TIME_COMM);
      }

      if (Lmask) boundary_clear_mask(iset);

      timer->stamp();

      // sectors but no colors (could also be no sectors)
      // random selection of sites in iset

      if (sweepflag == RANDOM) {
	site2i = set[iset].site2i;
	nrange = set[iset].nlocal;
	nselect = set[iset].nselect;
	for (i = 0; i < nselect; i++) 
	  sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
	(this->*sweep)(nselect,sitelist);
	nattempt += nselect;

      // sectors but no colors, or colors but no sectors
      // ordered sweep over all sites in iset

      } else if (bothflag == 0) {
	for (i = 0; i < set[iset].nloop; i++)
	  (this->*sweep)(set[iset].nlocal,set[iset].site2i);
	nattempt += set[iset].nselect;

      // sectors and colors
      // icolor loop is over all colors in a sector
      // jset = set that contains sites of one color in one sector
      // ordered sweep over all sites in jset

      } else {
	for (icolor = 0; icolor < ncolors; icolor++) {
	  jset = nsector + iset*ncolors + icolor;
	  for (i = 0; i < set[jset].nloop; i++)
	    (this->*sweep)(set[jset].nlocal,set[jset].site2i);
	  nattempt += set[jset].nselect;
	}
      }

      timer->stamp(TIME_SOLVE);

      if (nprocs > 1) {
	if (sectorflag) comm->reverse_sector(iset);
	else comm->all_reverse();
	timer->stamp(TIME_COMM);
      }
    }

    if (allow_update) user_update(dt_rkmc);

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ----------------------------------------------------------------------
   iterate the user_update routine only
   app is responsible for doing communciation in user_update()
------------------------------------------------------------------------- */

void AppLattice::iterate_update_only(double stoptime,double dt)
{
  int done = 0;
  while (!done) {
    if (allow_update) user_update(dt);
    
    time += dt;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_nomask_nostrict(int n, int *list)
{
  for (int m = 0; m < n; m++)
    site_event_rejection(list[m],ranapp);
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_mask_nostrict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    if (mask[i]) continue;
    site_event_rejection(i,ranapp);
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_nomask_strict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    ranstrict->seed = siteseeds[i];
    site_event_rejection(i,ranstrict);
    siteseeds[i] = ranstrict->seed;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_mask_strict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    if (mask[i]) continue;
    ranstrict->seed = siteseeds[i];
    site_event_rejection(i,ranstrict);
    siteseeds[i] = ranstrict->seed;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::input_app(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sector(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal sector command");

  nsector_user = 0;
  if (strcmp(arg[0],"yes") == 0) sectorflag = 1;
  else if (strcmp(arg[0],"no") == 0) sectorflag = 0;
  else {
    sectorflag = 1;
    nsector_user = atoi(arg[0]);
    if (nsector_user != 2 && nsector_user != 4 && nsector_user != 8)
      error->all(FLERR,"Illegal sector command");
  }

  nstop = 1.0;
  tstop = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      nstop = atof(arg[iarg+1]);
      if (nstop <= 0.0) error->all(FLERR,"Illegal sector command");
      tstop = 0.0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      tstop = atof(arg[iarg+1]);
      if (tstop <= 0.0) error->all(FLERR,"Illegal sector command");
      nstop = 0.0;
      iarg += 2;
    } else error->all(FLERR,"Illegal sector command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sweep(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal sweep command");
  if (strcmp(arg[0],"random") == 0) sweepflag = RANDOM;
  else if (strcmp(arg[0],"raster") == 0) sweepflag = RASTER;
  else if (strcmp(arg[0],"color") == 0) sweepflag = COLOR;
  else if (strcmp(arg[0],"color/strict") == 0) sweepflag = COLOR_STRICT;
  else if (strcmp(arg[0],"none") == 0) sweepflag = NOSWEEP;
  else error->all(FLERR,"Illegal sweep command");

  Lmask = false;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mask") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sweep command");
      if (strcmp(arg[iarg+1],"no") == 0) Lmask = false;
      else if (strcmp(arg[iarg+1],"yes") == 0) Lmask = true;
      else error->all(FLERR,"Illegal sweep command");
      iarg += 2;
    } else error->all(FLERR,"Illegal sweep command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_update_only(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal update_only command");
  if (strcmp(arg[0],"yes") == 0) update_only = 1;    
  else if (strcmp(arg[0],"no") == 0) update_only = 0;
  else error->all(FLERR,"Illegal update_only command");

  if (update_only && !allow_update)
    error->all(FLERR,"App does not permit user_update yes");
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice::stats(char *strtmp)
{
  char big[8],format[64];
  strcpy(big,BIGINT_FORMAT);

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  if (solve) {
    sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
    sprintf(strtmp,format,time,naccept_all,0,nsweeps);
  } else {
    bigint nattempt_all;
    MPI_Allreduce(&nattempt,&nattempt_all,1,MPI_SPK_BIGINT,MPI_SUM,world);
    sprintf(format,"%%10g %%10%s %%10%s %%10d",&big[1],&big[1]);
    sprintf(strtmp,format,time,naccept_all,nattempt_all-naccept_all,nsweeps);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nreject","Nsweeps");
}

/* ----------------------------------------------------------------------
   create a subset of owned sites
   insure all ptrs in Set data struct are allocated or NULL
   isector = 0 = all sites (no sector)
   isector > 1 = sites within a sector
   icolor = 0 = all sites (no color)
   icolor > 1 = sites of a certain color
 ------------------------------------------------------------------------- */

void AppLattice::create_set(int iset, int isector, int icolor, Solve *oldsolve)
{
  // sector boundaries

  double xmid = 0.5 * (domain->subxlo + domain->subxhi);
  double ymid = 0.5 * (domain->subylo + domain->subyhi);
  double zmid = 0.5 * (domain->subzlo + domain->subzhi);

  // count sites in subset

  int flag,iwhich,jwhich,kwhich,msector,mcolor;

  int delcolor = delevent + delpropensity;

  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 2) msector = iwhich + 1;
      else if (nsector == 4) msector = 2*jwhich + iwhich + 1;
      else msector = 4*kwhich + 2*jwhich + iwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (icolor > 0) {
      mcolor = domain->lattice->id2color(id[i],delcolor);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) n++;
  }

  set[iset].nlocal = n;

  // setup site2i for sites in set

  memory->create(set[iset].site2i,n,"app:site2i");

  n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 2) msector = iwhich + 1;
      else if (nsector == 4) msector = 2*jwhich + iwhich + 1;
      else msector = 4*kwhich + 2*jwhich + iwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (icolor > 0) {
      mcolor = domain->lattice->id2color(id[i],delcolor);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) set[iset].site2i[n++] = i;
  }

  // setup i2site for sites in set, only for KMC solver
  // i2site = 0 to nsite-1 for owned points in set, else -1

  if (solve) {
    (int *) memory->create(set[iset].i2site,nlocal+nghost,"app:i2site");
    for (int i = 0; i < nlocal+nghost; i++) set[iset].i2site[i] = -1;
    for (int i = 0; i < set[iset].nlocal; i++) 
      set[iset].i2site[set[iset].site2i[i]] = i;
  } else set[iset].i2site = NULL;

  // allocate propensity array for set

  memory->create(set[iset].propensity,n,"app:propensity");

  // allocate KMC solver for set
  // reuse old solver if provided, else delete it

  if (solve) {
    if (oldsolve) set[iset].solve = oldsolve;
    else set[iset].solve = solve->clone();
  } else {
    set[iset].solve = NULL;
    delete oldsolve;
  }

  // setup border arrays for set
  // nborder = # of sites in sector influenced by site outside sector
  // border = list of these sites, stored as lattice indices
  // bsites = scratch array for use by KMC solver
  // border is only used by KMC solver in sectors and masking
  // bsites is only used by KMC solver in sectors

  if ((solve && sectorflag) || Lmask) {
    set[iset].nborder = find_border_sites(iset);
    if (solve && sectorflag)
      memory->create(set[iset].bsites,set[iset].nborder,"app:bsites");
    else set[iset].bsites = NULL;
  } else {
    set[iset].nborder = 0;
    set[iset].border = NULL;
    set[iset].bsites = NULL;
  }
}

/* ----------------------------------------------------------------------
   free memory inside a set
   return Solver so caller can reuse it if desired
 ------------------------------------------------------------------------- */

Solve *AppLattice::free_set(int iset)
{
  memory->destroy(set[iset].border);
  memory->destroy(set[iset].bsites);
  memory->destroy(set[iset].propensity);
  memory->destroy(set[iset].site2i);
  memory->destroy(set[iset].i2site);
  return set[iset].solve;
}

/* ----------------------------------------------------------------------
   create list of border sites for a set
   border site = site in set with a 1 to Nlayer neighbor outside the set
   neighbor can be another owned site (outside set) or a ghost
   border = lattice index of the sites
 ------------------------------------------------------------------------- */

int AppLattice::find_border_sites(int isector)
{
  int i,j,m;

  int nlayer = delpropensity+delevent;
  int ntotal = nlocal + nghost;
  int nsites = set[isector].nlocal;
  int *site2i = set[isector].site2i;

  // flag sites with -1 that are not in sector
  // flag sites with 0 that are in sector

  int *flag;
  memory->create(flag,ntotal,"app:flag");
  for (i = 0; i < ntotal; i++) flag[i] = -1;
  for (m = 0; m < nsites; m++) flag[site2i[m]] = 0;

  // flag sector sites with -1 that have non-sector neighbor up to nlayer away

  for (int ilayer = 0; ilayer < nlayer; ilayer++) {
    for (m = 0; m < nsites; m++) {
      i = site2i[m];
      if (flag[i]) continue;
      for (j = 0; j < numneigh[i]; j++)
	if (flag[neighbor[i][j]] < 0) break;
      if (j < numneigh[i]) flag[i] = 1;
    }
    for (m = 0; m < nsites; m++) {
      i = site2i[m];
      if (flag[i] > 0) flag[i] = -1;
    }
  }

  // nborder = # of border sites
  // allocate border and fill with site indices

  int nborder = 0;
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    if (flag[i] < 0) nborder++;
  }

  int *border;
  memory->create(border,nborder,"app:border");

  nborder = 0;
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    if (flag[i] < 0) border[nborder++] = i;
  }

  memory->destroy(flag);

  set[isector].border = border;
  return nborder;
}
  
/* ----------------------------------------------------------------------
   unset all mask values of owned sites in iset whose propensity
     could change due to events on sites one neighbor outside the set
   border list stores indices of these sites
   their mask value may be out of date, due to state change in other sets
 ------------------------------------------------------------------------- */

void AppLattice::boundary_clear_mask(int iset)
{
  int *border = set[iset].border;
  int nborder = set[iset].nborder;

  for (int m = 0; m < nborder; m++) mask[border[m]] = 0;
}

/* ----------------------------------------------------------------------
   push new site onto stack and assign new id
 ------------------------------------------------------------------------- */

void AppLattice::push_new_site(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  // This can be used to screen out unwanted spin values
  // int isite = iarray[0][i];

  cluststack->push(i);
  cluster_ids[i] = id;
}

/* ----------------------------------------------------------------------
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */

void AppLattice::push_connected_neighbors(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  int ii;
  int isite = iarray[0][i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && cluster_ids[ii] == 0) {
      cluststack->push(ii);
      cluster_ids[ii] = id;
    }
  }
}

/* ----------------------------------------------------------------------
   add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */

void AppLattice::connected_ghosts(int i, int* cluster_ids, 
				  Cluster* clustlist, int idoffset)
{
  int iclust;
  int ii;
  int isite = iarray[0][i];

  // check if this was a site that was ignored

  if (cluster_ids[i] == 0) return;

  iclust = cluster_ids[i]-idoffset;
  
  // add ghost cluster to neighbors of local cluster

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && ii >= nlocal) {
      clustlist[iclust].add_neigh(cluster_ids[ii]);
    }
  }
}

/* ----------------------------------------------------------------------
   grow per-site arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AppLattice::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  memory->grow(id,nmax,"app:id");
  memory->grow(xyz,nmax,3,"app:xyz");
  memory->grow(owner,nmax,"app:owner");
  memory->grow(index,nmax,"app:index");

  memory->grow(numneigh,nmax,"app:numneigh");
  if (maxneigh) memory->grow(neighbor,nmax,maxneigh,"app:neighbor");

  for (int i = 0; i < ninteger; i++)
    memory->grow(iarray[i],nmax,"app:iarray");
  for (int i = 0; i < ndouble; i++)
    memory->grow(darray[i],nmax,"app:darray");

  grow_app();
}

/* ----------------------------------------------------------------------
   add an owned site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void AppLattice::add_site(tagint n, double x, double y, double z)
{
  if (nlocal == nmax) grow(0);

  id[nlocal] = n;
  xyz[nlocal][0] = x;
  xyz[nlocal][1] = y;
  xyz[nlocal][2] = z;

  owner[nlocal] = me;
  index[nlocal] = nlocal;

  for (int i = 0; i < ninteger; i++) iarray[i][nlocal] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][nlocal] = 0.0;

  nlocal++;
}

/* ----------------------------------------------------------------------
   add a ghost site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void AppLattice::add_ghost(tagint n, double x, double y, double z,
			   int proc, int index_owner)
{
  if (nlocal+nghost == nmax) grow(0);

  int m = nlocal + nghost;

  id[m] = n;
  xyz[m][0] = x;
  xyz[m][1] = y;
  xyz[m][2] = z;

  owner[m] = proc;
  index[m] = index_owner;

  for (int i = 0; i < ninteger; i++) iarray[i][m] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][m] = 0.0;

  nghost++;
}

/* ----------------------------------------------------------------------
   set neighbor connectivity for owned site I
   nvalues = # of neighbors
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppLattice::add_neighbors(int i, int nvalues, char **values)
{
  numneigh[i] = nvalues;
  for (int m = 0; m < nvalues; m++)
    neighbor[i][m] = atoi(values[m]);
}

/* ----------------------------------------------------------------------
   set values for owned site I
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppLattice::add_values(int i, char **values)
{
  for (int m = 0; m < ninteger; m++) iarray[m][i] = atoi(values[m]);
  for (int m = 0; m < ndouble; m++) darray[m][i] = atof(values[m+ninteger]);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   lo,hi = inclusive bounds
   5 possibilities:
     (1) i = i to i, (2) * = lo to hi,
     (3) i* = i to hi, (4) *j = lo to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void AppLattice::bounds(char *str, int lo, int hi, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(str),hi);
  } else if (strlen(str) == 1) {
    nlo = lo;
    nhi = hi;
  } else if (ptr == str) {
    nlo = lo;
    nhi = MIN(atoi(ptr+1),hi);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),lo);
    nhi = hi;
  } else {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(ptr+1),hi);
  }
}

/* ----------------------------------------------------------------------
   print connectivity stats
 ------------------------------------------------------------------------- */

void AppLattice::print_connectivity()
{
  int i;

  tagint min = maxneigh;
  tagint max = 0;

  for (i = 0; i < nlocal; i++) {
    min = MIN(min,numneigh[i]);
    max = MAX(max,numneigh[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_SPK_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_SPK_TAGINT,MPI_MAX,world);

  tagint *count = new tagint[maxall+1];
  tagint *countall = new tagint[maxall+1];

  for (i = 0; i <= maxall; i++) count[i] = 0;

  for (i = 0; i < nlocal; i++) count[numneigh[i]]++;
  MPI_Allreduce(count,countall,maxall+1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (me == 0)
    for (i = minall; i <= maxall; i++) {
      if (screen)
	fprintf(screen,"  " TAGINT_FORMAT " sites have %d neighbors\n",\
		countall[i],i);
      if (logfile)
	fprintf(logfile,"  " TAGINT_FORMAT " sites have %d neighbors\n",
		countall[i],i);
    }

  delete [] count;
  delete [] countall;
}
