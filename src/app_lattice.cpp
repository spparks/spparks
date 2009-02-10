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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "finish.h"
#include "cluster.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

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

  temperature = 0.0;
  
  lattice = NULL;
  iarray = NULL;
  darray = NULL;
  ninteger = ndouble = 0;

  propensity = NULL;
  i2site = NULL;

  comm = NULL;
  sweep = NULL;

  dt_sweep = 0.0;
  naccept = nattempt = 0;
  nsweeps = 0;

  // app can override these values in its constructor

  delpropensity = 1;
  delevent = 0;
  numrandom = 1;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  for (int i = 0; i < nset; i++) {
    memory->sfree(set[i].border);
    memory->sfree(set[i].bsites);
    delete set[i].solve;
    memory->sfree(set[i].propensity);
    memory->sfree(set[i].site2i);
    memory->sfree(set[i].i2site);
  }
  delete [] set;

  delete ranapp;
  delete ranstrict;
  memory->sfree(siteseeds);
  memory->sfree(sitelist);
  memory->sfree(mask);

  memory->sfree(lattice);
  for (int i = 0; i < ninteger; i++) memory->sfree(iarray[i]);
  for (int i = 0; i < ndouble; i++) memory->sfree(darray[i]);
  delete [] iarray;
  delete [] darray;

  delete comm;

  delete [] latfile;
  delete [] infile;

  memory->sfree(id);
  memory->sfree(owner);
  memory->sfree(index);
  memory->destroy_2d_T_array(xyz);

  memory->sfree(numneigh);
  memory->destroy_2d_T_array(neighbor);
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  // error checks

  if (solve == NULL && sweepflag == NOSWEEP)
    error->all("Lattice app needs a KMC or rejection KMC solver");
  if (solve && sweepflag != NOSWEEP)
    error->all("Lattice app cannot use both a KMC and rejection KMC solver");

  if (solve && allow_kmc == 0)
    error->all("KMC events are not implemented in app");
  if (sweepflag != NOSWEEP && allow_rejection == 0)
    error->all("Rejection events are not implemented in app");
  if (sweepflag != NOSWEEP && Lmask && allow_masking == 0)
    error->all("Mask logic not implemented in app");

  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all("Cannot use KMC solver in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
    error->all("Cannot use random rejection KMC in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
    error->all("Cannot use raster rejection KMC in parallel with no sectors");
  if (sectorflag && (sweepflag == COLOR || sweepflag == COLOR_STRICT))
    error->all("Cannot use colored rejection KMC with sectors");

  if (sweepflag && dt_sweep == 0.0)
    error->all("Lattice app did not set dt_sweep");

  // app-specific initialization

  init_app();

  // if sectors, determine number of sectors

  int nsector = 1;
  if (sectorflag) {
    if (dimension == 2) nsector = 4;
    else nsector = 8;
  }

  // if coloring, determine number of colors
  // setup test for create_set
  // check periodicity against lattice extent

  int ncolors = 1;
  if (sweepflag == COLOR || sweepflag == COLOR_STRICT) {
    int delcolor = delevent + delpropensity;
    if (latstyle == SQ_4N) {
      if (delcolor == 1) ncolors = 2;
      if (nx % 2 || ny % 2)
	error->all("Color stencil is incommensurate with lattice size");
    } else if (latstyle == SQ_8N) {
      ncolors = (delcolor+1)*(delcolor+1);
      if (nx % (delcolor+1) || ny % (delcolor+1))
	error->all("Color stencil is incommensurate with lattice size");
    } else if (latstyle == TRI) {
      if (delcolor == 1) ncolors = 4;
      if (nx % 2)
	error->all("Color stencil is incommensurate with lattice size");
    } else if (latstyle == SC_6N) {
      if (delcolor == 1) ncolors = 2;
      if (nx % 2 || ny % 2 || nz % 2)
	error->all("Color stencil is incommensurate with lattice size");
    } else if (latstyle == SC_26N) {
      ncolors = (delcolor+1)*(delcolor+1)*(delcolor+1);
      if (nx % (delcolor+1) || ny % (delcolor+1) || nz % (delcolor+1))
	error->all("Color stencil is incommensurate with lattice size");
    } else if (latstyle == FCC) {
      if (delcolor == 1) ncolors = 4;
    } else if (latstyle == BCC) {
      if (delcolor == 1) ncolors = 2;
    }

    if (ncolors == 1)
      error->all("Cannot color this combination of lattice and app");
  }

  // create sets based on sectors and coloring
  // only do this on first init

  if (set == NULL) {
    if (nsector == 1 && ncolors == 1) {
      nset = 1;
      set = new Set[nset];
      create_set(0,0,0);
    } else if (nsector > 1) {
      nset = nsector;
      set = new Set[nset];
      for (int i = 0; i < nset; i++) create_set(i,i+1,0);
    } else if (ncolors > 1) {
      nset = ncolors;
      set = new Set[nset];
      for (int i = 0; i < nset; i++) create_set(i,0,i+1);
    }
  }

  // initialize mask array

  if (!Lmask && mask) {
    memory->sfree(mask);
    mask = NULL;
  }
  if (Lmask && mask == NULL) {
    mask = (char *)
      memory->smalloc((nlocal+nghost)*sizeof(char),"applattice:mask");
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
    memory->sfree(siteseeds);
    ranstrict = NULL;
    siteseeds = NULL;
  }

  if (sweepflag == COLOR_STRICT && ranstrict == NULL) {
    ranstrict = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    siteseeds = 
      (int *) memory->smalloc(nlocal*sizeof(int),"applattice:siteseeds");
    for (int i = 0; i < nlocal; i++) {
      ranstrict->reset(seed,id[i],100);
      siteseeds[i] = ranstrict->seed;
    }
  }

  // initialize comm, both for this proc's full domain and sectors
  // only do this on first init

  if (comm == NULL) {
    comm = new CommLattice(spk);
    comm->init(nsector,delpropensity,delevent,NULL);
  }

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

  if (sweepflag != NOSWEEP) {
    if (sweepflag == RANDOM) {
      if (nstop > 0.0) {
	for (int i = 0; i < nset; i++)
	  set[i].nselect = static_cast<int> (nstop*set[i].nlocal);
      }
      if (tstop > 0.0) {
	double n = tstop / (dt_sweep/nglobal);
	for (int i = 0; i < nset; i++)
	  set[i].nselect = static_cast<int> (n/nglobal * set[i].nlocal);
      }

    } else if (sweepflag == RASTER || 
	       sweepflag == COLOR || sweepflag == COLOR_STRICT) {
      int n;
      if (nstop > 0.0) n = static_cast<int> (nstop);
      if (tstop > 0.0) n = static_cast<int> (tstop/dt_sweep);
      for (int i = 0; i < nset; i++) {
	set[i].nloop = n;
	set[i].nselect = n * set[i].nlocal;
      }
    }

    double nme = 0.0;
    for (int i = 0; i < nset; i++) nme += set[i].nselect;
    double ntotal;
    MPI_Allreduce(&nme,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);

    dt_rkmc = ntotal/nglobal * dt_sweep;
    if (dt_rkmc == 0.0)
      error->all("Choice of sector stop led to no rKMC events");
    dt_rkmc = MIN(dt_rkmc,stoptime-time);
  }


  // setup sitelist if sweepflag = RANDOM
  // do this every init since sector timestep could have changed

  if (sweepflag == RANDOM) {
    memory->sfree(sitelist);
    int n = 0;
    for (int i = 0; i < nset; i++) n = MAX(n,set[i].nselect);
    sitelist = (int *) memory->smalloc(n*sizeof(int),"applattice:sitelist");
  }

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

  // initialize output

  output->init(time);
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");

  // time_eps eliminates overruns due to finite machine precision
  // can be set to zero by an app that wants to ignore it

  timer->init();
  stoptime = time + atof(arg[0]) - time_eps;
  if (time >= stoptime) {
    Finish finish(spk);
    return;
  }

  // initialize and run

  init();

  timer->barrier_start(TIME_LOOP);

  if (solve) {
    if (sectorflag == 0) iterate_kmc_global(stoptime);
    else iterate_kmc_sector(stoptime);
  } else iterate_rejection(stoptime);

  timer->barrier_stop(TIME_LOOP);
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   KMC solver on entire domain
   can only be invoked in serial
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_global(double stoptime)
{
  int isite;
  double dt;
  
  // global KMC runs with one set
  // save ptr to system solver

  Solve *hold_solve = solve;
  solve = set[0].solve;
  propensity = set[0].propensity;
  i2site = set[0].i2site;

  int done = 0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    if (isite < 0) done = 1;
    else {
      naccept++;
      site_event(isite,ranapp);
      time += dt;
      timer->stamp(TIME_APP);
    }

    if (time >= stoptime) done = 1;
    output->compute(time,done);
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

	solve->update(0,NULL,NULL);

	timer->stamp(TIME_SOLVE);
	
	if (isite < 0) done = 1;
	else {
	  naccept++;
	  timesector += dt;	
	  if (timesector >= dt_kmc) done = 1;
	  else site_event(site2i[isite],ranapp);
	  timer->stamp(TIME_APP);
	}
      }
      
      if (nprocs > 1) {
	comm->reverse_sector(iset);
	timer->stamp(TIME_COMM);
      }
    }

    // keep looping until overall time threshhold reached
    
    time += dt_kmc;
    if (time >= stoptime) alldone = 1;

    output->compute(time,alldone);
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
  int i,nselect,nrange;
  int *site2i;

  int done = 0;
  while (!done) {
    for (int iset = 0; iset < nset; iset++) {
      if (nprocs > 1) {
	timer->stamp();
	if (sectorflag) comm->sector(iset);
	else comm->all();
	timer->stamp(TIME_COMM);
      }

      if (Lmask) boundary_clear_mask(iset);

      timer->stamp();

      if (sweepflag == RANDOM) {
	site2i = set[iset].site2i;
	nrange = set[iset].nlocal;
	nselect = set[iset].nselect;
	for (i = 0; i < nselect; i++) 
	  sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
	(this->*sweep)(nselect,sitelist);
	nattempt += nselect;

      } else {
	for (i = 0; i < set[iset].nloop; i++)
	  (this->*sweep)(set[iset].nlocal,set[iset].site2i);
	nattempt += set[iset].nselect;
      }

      timer->stamp(TIME_SOLVE);

      if (nprocs > 1) {
	if (sectorflag) comm->reverse_sector(iset);
	else comm->all_reverse();
	timer->stamp(TIME_COMM);
      }

    }

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;

    output->compute(time,done);
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

void AppLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"sweep") == 0) set_sweep(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) output->set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice::input_app(char *command, int narg, char **arg)
{
  error->all("Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sector(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal sector command");
  if (strcmp(arg[0],"yes") == 0) sectorflag = 1;
  else if (strcmp(arg[0],"no") == 0) sectorflag = 0;
  else error->all("Illegal sector command");

  nstop = 1.0;
  tstop = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nstop") == 0) {
      if (iarg+2 > narg) error->all("Illegal sector command");
      nstop = atof(arg[iarg+1]);
      if (nstop <= 0.0) error->all("Illegal sector command");
      tstop = 0.0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tstop") == 0) {
      if (iarg+2 > narg) error->all("Illegal sector command");
      tstop = atof(arg[iarg+1]);
      if (tstop <= 0.0) error->all("Illegal sector command");
      nstop = 0.0;
      iarg += 2;
    } else error->all("Illegal sector command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sweep(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal sweep command");
  if (strcmp(arg[0],"random") == 0) sweepflag = RANDOM;
  else if (strcmp(arg[0],"raster") == 0) sweepflag = RASTER;
  else if (strcmp(arg[0],"color") == 0) sweepflag = COLOR;
  else if (strcmp(arg[0],"color/strict") == 0) sweepflag = COLOR_STRICT;
  else error->all("Illegal sweep command");

  Lmask = false;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mask") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep command");
      if (arg[iarg+1],"no" == 0) Lmask = false;
      else if (arg[iarg+1],"yes" == 0) Lmask = true;
      else error->all("Illegal sweep command");
      iarg += 2;
    } else error->all("Illegal sweep command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice::stats(char *strtmp)
{
  int naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_INT,MPI_SUM,world);
  if (solve) sprintf(strtmp,"%10g %10d %10d %10d",time,naccept_all,0,0);
  else {
    int nattempt_all;
    MPI_Allreduce(&nattempt,&nattempt_all,1,MPI_INT,MPI_SUM,world);
    sprintf(strtmp,"%10g %10d %10d %10d",
	    time,naccept_all,nattempt_all-naccept_all,nsweeps);
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

void AppLattice::create_set(int iset, int isector, int icolor)
{
  // sector boundaries

  double xmid = 0.5 * (subxlo + subxhi);
  double ymid = 0.5 * (subylo + subyhi);
  double zmid = 0.5 * (subzlo + subzhi);

  // count sites in subset

  int flag,iwhich,jwhich,kwhich,msector,mcolor;

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

      if (dimension == 2) msector = 2*iwhich + jwhich + 1;
      else msector = 4*iwhich + 2*jwhich + kwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (icolor > 0) {
      mcolor = id2color(id[i]);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) n++;
  }

  set[iset].nlocal = n;

  // setup site2i for sites in set

  set[iset].site2i =
    (int *) memory->smalloc(n*sizeof(int),"applattice:site2i");

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

      if (dimension == 2) msector = 2*iwhich + jwhich + 1;
      else msector = 4*iwhich + 2*jwhich + kwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (icolor > 0) {
      mcolor = id2color(id[i]);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) set[iset].site2i[n++] = i;
  }

  // setup i2site for sites in set, only for KMC solver
  // i2site = 0 to nsite-1 for owned points in set, else -1

  if (solve) {
    set[iset].i2site =
      (int *) memory->smalloc((nlocal+nghost)*sizeof(int),"applattice:i2site");
      for (int i = 0; i < nlocal+nghost; i++) set[iset].i2site[i] = -1;
      for (int i = 0; i < set[iset].nlocal; i++) 
	set[iset].i2site[set[iset].site2i[i]] = i;
  } else set[iset].i2site = NULL;

  // allocate propensity array for set

  set[iset].propensity =
    (double *) memory->smalloc(n*sizeof(double),"applattice:propensity");

  // allocate KMC solver for set

  if (solve) set[iset].solve = solve->clone();
  else set[iset].solve = NULL;

  // setup border arrays for set
  // nborder = # of sites in sector influenced by site outside sector
  // border = list of these sites, stored as lattice indices
  // bsites = scratch array for use by KMC solver
  // border is only used by KMC solver in sectors and masking
  // bsites is only used by KMC solver in sectors

  if ((solve && sectorflag) || Lmask) {
    set[iset].nborder = find_border_sites(iset);
    if (solve && sectorflag)
      set[iset].bsites =
	(int *) memory->smalloc(set[iset].nborder*sizeof(int),
				"applattice:bsites");
    else set[iset].bsites = NULL;
  } else {
    set[iset].nborder = 0;
    set[iset].border = NULL;
    set[iset].bsites = NULL;
  }
}

/* ----------------------------------------------------------------------
   convert a lattice ID (1 to Nsites) to a color (1 to Ncolor)
------------------------------------------------------------------------- */

int AppLattice::id2color(int idsite)
{
  int i,j,k,ncolors,ncolor1d,icolor;

  idsite--;
  int delcolor = delevent + delpropensity;

  if (latstyle == SQ_4N) {
    i = idsite % nx;
    j = idsite / nx;
    icolor = (i+j) % 2;

  } else if (latstyle == SQ_8N) {
    ncolor1d = delcolor+1;
    i = idsite % nx;
    j = idsite / nx;
    icolor = ncolor1d*(j%ncolor1d) + i%ncolor1d;

  } else if (latstyle == TRI) {
    icolor = idsite % 4;

  } else if (latstyle == SC_6N) {
    i = idsite % nx;
    j = (idsite%(nx*ny)) / nx;
    k = idsite / (nx*ny);
    icolor = (i+j+k) % 2;

  } else if (latstyle == SC_26N) {
    ncolor1d = delcolor+1;
    i = idsite % nx;
    j = (idsite%(nx*ny)) / nx;
    k = idsite / (nx*ny);
    icolor = ncolor1d*ncolor1d*(k%ncolor1d) + 
      ncolor1d*(j%ncolor1d) + i%ncolor1d;

  } else if (latstyle == FCC) {
    icolor = idsite % 4;

  } else if (latstyle == BCC) {
    icolor = idsite % 2;
  }

  return icolor+1;
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

  int *flag = (int *) memory->smalloc(ntotal*sizeof(int),"applattice:flag");
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

  int *border = (int *)
    memory->smalloc(nborder*sizeof(int),"applattice:border");

  nborder = 0;
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    if (flag[i] < 0) border[nborder++] = i;
  }

  memory->sfree(flag);

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
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */

void AppLattice::push_connected_neighbors(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  int ii;
  int isite = lattice[i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (lattice[ii] == isite && cluster_ids[ii] == 0) {
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
  int isite = lattice[i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (lattice[ii] == isite && ii >= nlocal) {
      iclust = cluster_ids[i]-idoffset;
      // Add ghost cluster to neighbors of local cluster
      clustlist[iclust].add_neigh(cluster_ids[ii]);
    }
  }
}
