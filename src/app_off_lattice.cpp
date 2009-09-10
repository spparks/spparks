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
#include "string.h"
#include "stdlib.h"
#include "app_off_lattice.h"
#include "comm_off_lattice.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 100
#define MAXNEIGH 1000

enum{NOSWEEP,RANDOM,RASTER};
enum{INTERIOR,EDGE,GHOST};            // same as in comm_off_lattice.cpp

/* ---------------------------------------------------------------------- */

AppOffLattice::AppOffLattice(SPPARKS *spk, int narg, char **arg) : 
  App(spk,narg,arg)
{
  appclass = OFF_LATTICE;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  sectorflag = 0;
  nset = 0;
  set = NULL;
  nstop = 1.0;
  tstop = 0.0;

  sweepflag = NOSWEEP;
  ranapp = NULL;
  sitelist = NULL;

  temperature = 0.0;
  
  site = NULL;
  iarray = NULL;
  darray = NULL;
  ninteger = ndouble = 0;

  nlocal = nghost = nmax = 0;
  nfree = 0;
  freehead = -1;

  id = NULL;
  xyz = NULL;
  bin = NULL;
  next = prev = nextimage = NULL;

  comm = NULL;
  //sweep = NULL;

  binhead = binflag = NULL;
  binoffset = NULL;
  nbinimages = NULL;
  binimage = NULL;

  nstencil = 0;
  stencil = NULL;

  neighs = new int[MAXNEIGH];

  naccept = nattempt = 0;
  nsweeps = 0;
}

/* ---------------------------------------------------------------------- */

AppOffLattice::~AppOffLattice()
{
  for (int i = 0; i < nset; i++) {
    delete set[i].solve;
    memory->sfree(set[i].propensity);
    memory->sfree(set[i].site2i);
    memory->sfree(set[i].i2site);
  }
  delete [] set;

  delete ranapp;
  memory->sfree(sitelist);

  memory->sfree(site);
  for (int i = 0; i < ninteger; i++) memory->sfree(iarray[i]);
  for (int i = 0; i < ndouble; i++) memory->sfree(darray[i]);
  delete [] iarray;
  delete [] darray;

  memory->sfree(id);
  memory->destroy_2d_T_array(xyz);
  memory->sfree(bin);
  memory->sfree(next);
  memory->sfree(prev);
  memory->sfree(nextimage);

  delete comm;

  memory->sfree(binhead);
  memory->sfree(binflag);
  memory->destroy_2d_int_array(binoffset);
  memory->sfree(nbinimages);
  memory->destroy_2d_int_array(binimage);

  delete [] stencil;
  delete [] neighs;
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"sweep") == 0) set_sweep(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) output->add_dump(narg,arg);
  else if (strcmp(command,"dump_one") == 0) output->dump_one(narg,arg,time);
  else if (strcmp(command,"dump_modify") == 0) output->dump_modify(narg,arg);
  else if (strcmp(command,"undump") == 0) output->undump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::init()
{
  int i,j,k,m,ix,iy,iz;

  // error checks

  if (solve == NULL && sweepflag == NOSWEEP)
    error->all("App needs a KMC or rejection KMC solver");
  if (solve && sweepflag != NOSWEEP)
    error->all("App cannot use both a KMC and rejection KMC solver");

  if (solve && allow_kmc == 0)
    error->all("KMC events are not implemented in app");
  if (sweepflag != NOSWEEP && allow_rejection == 0)
    error->all("Rejection events are not implemented in app");

  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all("Cannot use KMC solver in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
    error->all("Cannot use random rejection KMC in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
    error->all("Cannot use raster rejection KMC in parallel with no sectors");

  if (sweepflag && dt_sweep == 0.0) error->all("App did not set dt_sweep");

  // if sectors, set number of sectors

  nsector = 1;
  if (sectorflag) {
    if (nsector_user) nsector = nsector_user;
    else if (dimension == 1) nsector = 2;
    else if (dimension == 2) nsector = 4;
    else nsector = 8;

    if (dimension == 3) {
      if (nsector == 2 && (ny_procs != 1 || nz_procs != 1))
	error->all("Invalid number of sectors");
      if (nsector == 4 && nz_procs != 1)
	error->all("Invalid number of sectors");
    }
    if (dimension == 2) {
      if (nsector == 2 && ny_procs != 1)
	error->all("Invalid number of sectors");
      if (nsector == 8)
	error->all("Invalid number of sectors");
    }
    if (dimension == 1 && nsector != 2)
      error->all("Invalid number of sectors");
  }

  // create sets based on sectors
  // only do this on first init

  if (set == NULL) {
    if (nsector == 1) {
      nset = 1;
      set = new Set[nset];
      create_set(0,0);
    } else if (nsector > 1) {
      nset = nsector;
      set = new Set[nset];
      for (int i = 0; i < nset; i++) create_set(i,i+1);
    }
  }

  // setup ranapp RN generator, only on first init
  // setup ranapp so different on every proc

  if (ranapp == NULL) {
    ranapp = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    ranapp->reset(seed,me,100);
  }

  // app-specific initialization, after general initialization

  init_app();

  // create bins for sites
  // do this on every init in case cutoff changed
  // bin size must be >= cutoff set by parameters in init_app()
  // require 2 bins within sub-domain if sectoring in a dimension
  // else require 1 bin within sub-domain
  // require an even number of bins if sectoring in a dimension
  // add 2 bins for ghost regions

  double cutoff = delpropensity + delevent;
  nbinx = static_cast<int> (xprd/nx_procs/cutoff);
  nbiny = nbinz = 1;
  if (dimension >= 2) nbiny = static_cast<int> (yprd/ny_procs/cutoff);
  if (dimension == 3) nbinz = static_cast<int> (zprd/nz_procs/cutoff);

  if (nsector >= 2) nbinx = nbinx/2 * 2;
  if (nsector >= 4) nbiny = nbiny/2 * 2;
  if (nsector == 8) nbinz = nbinz/2 * 2;

  if (nsector >= 2 && nbinx < 2)
    error->all("Application cutoff is too big for processor sub-domain");
  if (nsector >= 4 && nbiny < 2)
    error->all("Application cutoff is too big for processor sub-domain");
  if (nsector == 8 && nbinz < 2)
    error->all("Application cutoff is too big for processor sub-domain");
  if (nbinx < 1 || nbiny < 1 || nbinz < 1)
    error->all("Application cutoff is too big for processor sub-domain");

  binx = xprd / nbinx;
  biny = yprd / nbiny;
  binz = zprd / nbinz;
  invbinx = 1.0 / binx;
  invbiny = 1.0 / biny;
  invbinz = 1.0 / binz;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  Bins: %d %d %d\n",nbinx,nbiny,nbinz);
      fprintf(screen,"  Bin sizes: %g %g %g\n",binx,biny,binz);
    }
    if (logfile) {
      fprintf(logfile,"  Bins: %d %d %d\n",nbinx,nbiny,nbinz);
      fprintf(logfile,"  Bin sizes: %g %g %g\n",binx,biny,binz);
    }
  }

  nbinx += 2;
  if (dimension >= 2) nbiny += 2;
  if (dimension == 3) nbinz += 2;

  nbins = nbinx*nbiny*nbinz;

  memory->sfree(binhead);
  memory->sfree(binflag);
  memory->destroy_2d_int_array(binoffset);
  memory->sfree(nbinimages);
  memory->destroy_2d_int_array(binimage);

  binhead = (int *) memory->smalloc(nbins*sizeof(int),"app:binhead");
  binflag = (int *) memory->smalloc(nbins*sizeof(int),"app:binflag");
  binoffset = memory->create_2d_int_array(nbins,3,"app:binoffset");
  nbinimages = (int *) memory->smalloc(nbins*sizeof(int),"app:nbinimages");
  if (dimension == 3)
    binimage = memory->create_2d_int_array(nbins,7,"app:binimage");
  else if (dimension == 2)
    binimage = memory->create_2d_int_array(nbins,3,"app:binimage");
  else
    binimage = memory->create_2d_int_array(nbins,1,"app:binimage");

  if (dimension == 3) {
    for (i = 0; i < nbinx; i++) 
      for (j = 0; j < nbiny; j++) 
	for (k = 0; k < nbinz; k++) {
	  m = k*nbiny*nbinx + j*nbinx + i;
	  if (i == 0 || i == nbinx-1 || j == 0 || j == nbiny-1 || 
	      k == 0 || k == nbinz-1) binflag[m] = GHOST;
	  else if (i > 1 && i < nbinx-2 && j > 1 && j < nbiny-2 &&
		   k > 1 && k < nbinz-2) binflag[m] = INTERIOR;
	  else binflag[m] = EDGE;
	}
  } else if (dimension == 2) {
    for (i = 0; i < nbinx; i++) 
      for (j = 0; j < nbiny; j++) {
	m = j*nbinx + i;
	if (i == 0 || i == nbinx-1 || j == 0 || j == nbiny-1)
	  binflag[m] = GHOST;
	else if (i > 1 && i < nbinx-2 && j > 1 && j < nbiny-2)
	  binflag[m] = INTERIOR;
	else binflag[m] = EDGE;
      }
  } else {
    for (i = 0; i < nbinx; i++) 
      for (j = 0; j < nbiny; j++) {
	m = j*nbinx + i;
	if (i == 0 || i == nbinx-1 || j == 0 || j == nbiny-1)
	  binflag[m] = GHOST;
	else if (i > 1 && i < nbinx-2 && j > 1 && j < nbiny-2)
	  binflag[m] = INTERIOR;
	else binflag[m] = EDGE;
      }
  }

  for (m = 0; m < nbins; m++) {
    binhead[m] = -1;
    if (binflag[m] == GHOST) {
      ix = m % nbinx;
      iy = (m/nbinx) % nbiny;
      iz = m / (nbinx*nbiny);
      if (ix == 0) binoffset[m][0] = -1;
      else if (ix == nbinx-1) binoffset[m][0] = 1;
      else binoffset[m][0] = 0;
      if (dimension >= 2) {
	if (iy == 0) binoffset[m][1] = -1;
	else if (iy == nbiny-1) binoffset[m][1] = 1;
	else binoffset[m][1] = 0;
      } else binoffset[m][1] = 0;
      if (dimension == 3) {
	if (iz == 0) binoffset[m][2] = -1;
	else if (iz == nbinz-1) binoffset[m][2] = 1;
	else binoffset[m][2] = 0;
      } else binoffset[m][2] = 0;

    } else if (binflag[m] == EDGE) {
      nbinimages[m] = 0;
      ix = m % nbinx;
      iy = (m/nbinx) % nbiny;
      iz = m / (nbinx*nbiny);
      add_image_bins(m,ix,iy,iz);
    }
  }

  // create stencil for neighbor finding
  // do this on every init in case cutoff and bins changed

  delete [] stencil;
  if (dimension == 3) nstencil = 27;
  else nstencil = 9;
  stencil = new int[nstencil];
  nstencil = 0;
  if (dimension == 3) {
    for (k = -1; k <= 1; k++)
      for (j = -1; j <= 1; j++)
	for (i = -1; i <= 1; i++)
	  stencil[nstencil++] = k*nbiny*nbinx + j*nbinx + i;
  } else if (dimension == 2) {
    for (j = -1; j <= 1; j++)
      for (i = -1; i <= 1; i++)
	stencil[nstencil++] = j*nbinx + i;
  } else {
    for (i = -1; i <= 1; i++)
      stencil[nstencil++] = i;
  }

  // bin owned sites
  // loop in reverse order so linked list will be in forward order

  for (int i = nlocal-1; i >= 0; i--) {
    bin[i] = site2bin(i);
    add_to_bin(i,bin[i]);
  }

  // initialize comm, both for this proc's full domain and sectors
  // redo every run in case cutoff and bins changed

  delete comm;
  comm = new CommOffLattice(spk);
  comm->init(nsector);

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::setup()
{
  // app-specific setup, before propensities are computed

  setup_app();

  // comm ghost sites if needed

  if (nprocs == 1 && sectorflag == 0) comm->all();

  // convert rejection info to rKMC params
  // nloop and nselect are set whether sectoring is used or not

  if (sweepflag != NOSWEEP) {
    if (sweepflag == RANDOM) {
      if (nstop > 0.0) {
	for (int i = 0; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (nstop*set[i].nlocal);
	}
      }
      if (tstop > 0.0) {
	double n = tstop / (dt_sweep/nglobal);
	for (int i = 0; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (n/nglobal * set[i].nlocal);
	}
      }

    } else if (sweepflag == RASTER) {
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
  // do this every run since sector timestep could have changed

  if (sweepflag == RANDOM) {
    memory->sfree(sitelist);
    int n = 0;
    for (int i = 0; i < nset; i++) n = MAX(n,set[i].nselect);
    sitelist = (int *) memory->smalloc(n*sizeof(int),"app:sitelist");
  }

  // setup future output

  nextoutput = output->setup(time);
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::iterate()
{
  timer->barrier_start(TIME_LOOP);

  if (solve) {
    if (sectorflag == 0) iterate_kmc_global(stoptime);
    else iterate_kmc_sector(stoptime);
  } else iterate_rejection(stoptime);

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   KMC solver on entire domain
   can only be invoked in serial
 ------------------------------------------------------------------------- */

void AppOffLattice::iterate_kmc_global(double stoptime)
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

      // what does site_event do ??

      //site_event(isite);
      naccept++;
      time += dt;

      // what does update_event do ??

      //update_events();
      timer->stamp(TIME_APP);
    }

    if (time >= stoptime) done = 1;
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

void AppOffLattice::iterate_kmc_sector(double stoptime)
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

      //solve->update(nsites,bsites,propensity);
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
	    //site_event(site2i[isite]);
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

    // keep looping until overall time threshhold reached
    
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

void AppOffLattice::iterate_rejection(double stoptime)
{
  int i,nselect,nrange;
  int *site2i;

  int done = 0;
  while (!done) {
    for (int iset = 0; iset < nset; iset++) {
      /*
      if (nprocs > 1) {
	timer->stamp();
	if (bothflag) comm->sector(iset/ncolors);
	else if (sectorflag) comm->sector(iset);
	else comm->all();
	timer->stamp(TIME_COMM);
      }
      */

      timer->stamp();

      // sweep over random selection of sites in iset

      if (sweepflag == RANDOM) {
	site2i = set[iset].site2i;
	nrange = set[iset].nlocal;
	nselect = set[iset].nselect;
	for (i = 0; i < nselect; i++) 
	  sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
	for (int m = 0; m < nselect; m++)
	  site_event_rejection(sitelist[m],ranapp);

	nattempt += nselect;

      // ordered sweep

      } else {
	for (i = 0; i < set[iset].nloop; i++)
	  for (int m = 0; m < set[iset].nlocal; m++)
	    site_event_rejection(set[iset].site2i[m],ranapp);
	nattempt += set[iset].nselect;
      }

      timer->stamp(TIME_SOLVE);

      /*
      if (nprocs > 1) {
	if (bothflag) comm->reverse_sector(iset/ncolors);
	else if (sectorflag) comm->reverse_sector(iset);
	else comm->all_reverse();
	timer->stamp(TIME_COMM);
      }
      */
    }

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::input_app(char *command, int narg, char **arg)
{
  error->all("Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::set_sector(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal sector command");

  nsector_user = 0;
  if (strcmp(arg[0],"yes") == 0) sectorflag = 1;
  else if (strcmp(arg[0],"no") == 0) sectorflag = 0;
  else {
    sectorflag = 1;
    nsector_user = atoi(arg[0]);
    if (nsector_user != 2 && nsector_user != 4 && nsector_user != 8)
      error->all("Illegal sector command");
  }

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

void AppOffLattice::set_sweep(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal sweep command");
  if (strcmp(arg[0],"random") == 0) sweepflag = RANDOM;
  else if (strcmp(arg[0],"raster") == 0) sweepflag = RASTER;
  else error->all("Illegal sweep command");
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppOffLattice::stats(char *strtmp)
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

void AppOffLattice::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nreject","Nsweeps");
}

/* ----------------------------------------------------------------------
   create a subset of owned sites
   insure all ptrs in Set data struct are allocated or NULL
   isector = 0 = all sites (no sector)
   isector > 1 = sites within a sector
 ------------------------------------------------------------------------- */

void AppOffLattice::create_set(int iset, int isector)
{
  // sector boundaries

  double xmid = 0.5 * (subxlo + subxhi);
  double ymid = 0.5 * (subylo + subyhi);
  double zmid = 0.5 * (subzlo + subzhi);

  // NOTE: assume single proc for now

  int flag,iwhich,jwhich,kwhich,msector;

  int n = nlocal;
  set[iset].nlocal = n;

  // setup site2i for sites in set

  set[iset].site2i =
    (int *) memory->smalloc(n*sizeof(int),"app:site2i");

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
      else if (nsector == 4) msector = 2*iwhich + jwhich + 1;
      else msector = 4*iwhich + 2*jwhich + kwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (flag) set[iset].site2i[n++] = i;
  }

  // setup i2site for sites in set, only for KMC solver
  // i2site = 0 to nsite-1 for owned points in set, else -1

  if (solve) {
    set[iset].i2site =
      (int *) memory->smalloc((nlocal+nghost)*sizeof(int),"app:i2site");
      for (int i = 0; i < nlocal+nghost; i++) set[iset].i2site[i] = -1;
      for (int i = 0; i < set[iset].nlocal; i++) 
	set[iset].i2site[set[iset].site2i[i]] = i;
  } else set[iset].i2site = NULL;

  // allocate propensity array for set

  set[iset].propensity =
    (double *) memory->smalloc(n*sizeof(double),"app:propensity");

  // allocate KMC solver for set

  if (solve) set[iset].solve = solve->clone();
  else set[iset].solve = NULL;
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   name = "nlocal" or "lattice" or "xyz" or "iarrayN" or "darrayN"
   N in iarray and darray is an integer from 1 to ninteger or ndouble
 ------------------------------------------------------------------------- */

void *AppOffLattice::extract(char *name)
{
  if (strcmp(name,"nlocal") == 0) return (void *) &nlocal;
  if (strcmp(name,"xyz") == 0) return (void *) xyz;
  if (strstr(name,"iarray") == name) {
    int n = atoi(&name[6]);
    if (n < 1 || n > ninteger) return NULL;
    return (void *) iarray[n-1];
  }
  if (strstr(name,"darray") == name) {
    int n = atoi(&name[6]);
    if (n < 1 || n > ndouble) return NULL;
    return (void *) darray[n-1];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create list of neighbors of site I within distance cutoff
 ------------------------------------------------------------------------- */

void AppOffLattice::neighbor(int i, double cut)
{
  int j;
  double delx,dely,delz,rsq;

  double cutsq = cut*cut;

  double xtmp = xyz[i][0];
  double ytmp = xyz[i][1];
  double ztmp = xyz[i][2];

  numneigh = 0;

  int ibin = bin[i];

  for (int k = 0; k < nstencil; k++) {
    for (j = binhead[ibin+stencil[k]]; j >= 0; j = next[j]) {
      delx = xtmp - xyz[j][0];
      dely = ytmp - xyz[j][1];
      delz = ztmp - xyz[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (numneigh == MAXNEIGH) {   // delete this test eventually
	printf("NEIGH %d %d %d\n",i,ibin,numneigh);
	error->one("Too many neighbors per site");
      }
      if (rsq <= cutsq && i != j) neighs[numneigh++] = j;
    }
  }

  if (numneigh > MAXNEIGH) error->one("Too many neighbors per site");
}

/* ----------------------------------------------------------------------
   move site I to new position currently stored in xyz
   update all affected data structures
 ------------------------------------------------------------------------- */

void AppOffLattice::move(int i)
{
  int j;

  int oldbin = bin[i];
  int newbin = site2bin(i);

  // site stays in same bin
  // update its ghost image coords if necessary

  if (newbin == oldbin) {
    j = nextimage[i];
    while (j >= 0) {
      oldbin = bin[j];
      xyz[j][0] = xyz[i][0] + binoffset[oldbin][0]*xprd;
      xyz[j][1] = xyz[i][1] + binoffset[oldbin][1]*yprd;
      xyz[j][2] = xyz[i][2] + binoffset[oldbin][2]*zprd;
      j = nextimage[j];
    }
    return;
  }

  // site stays in INTERIOR
  // delete it from old bin, add it to new bin

  if (binflag[oldbin] == INTERIOR && binflag[newbin] == INTERIOR) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    return;
  }

  // site stays in EDGE
  // delete it from old bin, add it to new bin
  // delete ghost images
  // create new ghost images

  if (binflag[oldbin] == EDGE && binflag[newbin] == EDGE) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    delete_images(i);
    add_images(i,newbin);
    return;
  }

  // site moves from INTERIOR to EDGE
  // delete it from old bin, add it to new bin
  // create new ghost images

  if (binflag[oldbin] == INTERIOR && binflag[newbin] == EDGE) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    add_images(i,newbin);
    return;
  }

  // site moves from EDGE to INTERIOR
  // delete it from old bin, add it to new bin
  // delete ghost images

  if (binflag[oldbin] == EDGE && binflag[newbin] == INTERIOR) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    delete_images(i);
    return;
  }

  // site moves from EDGE to GHOST
  // map it back into simulation box via PBC, then is same as EDGE to EDGE
  // delete it from old bin, add it to new bin
  // delete ghost images
  // create new ghost images

  if (binflag[oldbin] == EDGE && binflag[newbin] == GHOST) {
    xyz[i][0] -= binoffset[newbin][0]*xprd;
    xyz[i][1] -= binoffset[newbin][1]*yprd;
    xyz[i][2] -= binoffset[newbin][2]*zprd;
    newbin = site2bin(i);
    if (binflag[newbin] != EDGE) error->one("PBC remap of site failed");
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    delete_images(i);
    add_images(i,newbin);
    return;
  }
}

/* ----------------------------------------------------------------------
   return bin ID containing site I
 ------------------------------------------------------------------------- */

int AppOffLattice::site2bin(int i)
{
  int ix,iy,iz;

  if (xyz[i][0] >= subxhi) ix = nbinx-1;
  else if (xyz[i][0] >= subxlo) 
    ix = static_cast<int> ((xyz[i][0]-subxlo)*invbinx) + 1;
  else ix = 0;

  if (dimension >= 2) {
    if (xyz[i][1] >= subyhi) iy = nbiny-1;
    else if (xyz[i][1] >= subylo) 
      iy = static_cast<int> ((xyz[i][1]-subylo)*invbiny) + 1;
    else iy = 0;
  } else iy = 0;

  if (dimension == 3) {
    if (xyz[i][2] >= subzhi) iz = nbinz-1;
    else if (xyz[i][2] >= subzlo) 
      iz = static_cast<int> ((xyz[i][2]-subzlo)*invbinz) + 1;
    else iz = 0;
  } else iz = 0;

  return iz*nbinx*nbiny + iy*nbinx + ix;
}

/* ----------------------------------------------------------------------
   delete site I from bin IBIN
 ------------------------------------------------------------------------- */

void AppOffLattice::delete_from_bin(int i, int ibin)
{
  if (prev[i] < 0) binhead[ibin] = next[i];
  else next[prev[i]] = next[i];
  if (next[i] >= 0) prev[next[i]] = prev[i];
}

/* ----------------------------------------------------------------------
   add site I to bin IBIN
 ------------------------------------------------------------------------- */

void AppOffLattice::add_to_bin(int i, int ibin)
{
  if (binhead[ibin] < 0) {
    binhead[ibin] = i;
    prev[i] = next[i] = -1;
  } else {
    prev[binhead[ibin]] = i;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
    prev[i] = -1;
  }
}

/* ----------------------------------------------------------------------
   delete ghost sites which are periodic images of site I
 ------------------------------------------------------------------------- */

void AppOffLattice::delete_images(int i)
{
  int oldbin;

  int j = nextimage[i];
  while (j >= 0) {
    oldbin = bin[j];
    delete_from_bin(j,oldbin);
    add_to_free(j);
    j = nextimage[j];
  }

  // only need to reset ptr of owned site

  nextimage[i] = -1;
}

/* ----------------------------------------------------------------------
   add ghost sites which are periodic images of site I in bin IBIN
 ------------------------------------------------------------------------- */

void AppOffLattice::add_images(int i, int ibin)
{
  int j,jbin;

  nextimage[i] = -1;
  int n = nbinimages[ibin];
  for (int m = 0; m < n; m++) {
    jbin = binimage[ibin][m];
    j = new_ghost();
    xyz[j][0] = xyz[i][0] + binoffset[jbin][0]*xprd;
    xyz[j][1] = xyz[i][1] + binoffset[jbin][1]*yprd;
    xyz[j][2] = xyz[i][2] + binoffset[jbin][2]*zprd;
    id[j] = id[i];
    bin[j] = jbin;
    add_to_bin(j,jbin);
    nextimage[j] = nextimage[i];
    nextimage[i] = j;
  }
}

/* ----------------------------------------------------------------------
   add ghost site I to free list
 ------------------------------------------------------------------------- */

void AppOffLattice::add_to_free(int i)
{
  next[i] = freehead;
  freehead = i;
  nfree++;
  nghost--;
}

/* ----------------------------------------------------------------------
   delete all ghost sites
   add them to free list
   set binhead of ghost bins to -1
 ------------------------------------------------------------------------- */

void AppOffLattice::delete_all_ghosts()
{
  if (nmax > nlocal) freehead = nlocal;
  for (int i = nlocal; i < nmax; i++) next[i] = i+1;
  if (nmax > nlocal) next[nmax-1] = -1;
  nfree = nmax-nlocal;
  nghost = 0;

  for (int i = 0; i < nbins; i++)
    if (binflag[i] == GHOST) binhead[i] = -1;
}

/* ----------------------------------------------------------------------
   return index of a new ghost site from free list
   grow per-site arrays if free list is empty
 ------------------------------------------------------------------------- */

int AppOffLattice::new_ghost()
{
  if (nfree == 0) {
    int oldmax = nmax;
    grow(0);
    add_free(oldmax);
  }

  int index = freehead;
  freehead = next[freehead];
  nfree--;
  nghost++;
  return index;
}

/* ----------------------------------------------------------------------
   add newly allocated sites to end of free list
   new sites are from oldmax to nmax
   assume either there were no free sites before allocation
   or there are no ghosts so that old free list ends at oldmax
 ------------------------------------------------------------------------- */

void AppOffLattice::add_free(int oldmax)
{
  if (freehead < 0) freehead = oldmax;
  else next[oldmax-1] = oldmax;
  for (int i = oldmax; i < nmax; i++) next[i] = i+1;
  next[nmax-1] = -1;
  nfree += nmax-oldmax;
}

/* ----------------------------------------------------------------------
   add periodic images of bin I,J,K to bin M
   calls add_images() recursively for image bins
   only add image bin if not already in binimage[] list
   can add up to 3 images in 2d, 7 in 3d
------------------------------------------------------------------------- */

void AppOffLattice::add_image_bins(int m, int i, int j, int k)
{
  int n;

  if (i == 1) {
    int ibin = k*nbinx*nbiny + j*nbinx + nbinx-1;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,nbinx-1,j,k);
  }
  if (i == nbinx-2) {
    int ibin = k*nbinx*nbiny + j*nbinx;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,0,j,k);
  }
  if (j == 1) {
    int ibin = k*nbinx*nbiny + (nbiny-1)*nbinx + i;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,i,nbiny-1,k);
  }
  if (j == nbiny-2) {
    int ibin = k*nbinx*nbiny + i;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,i,0,k);
  }
  if (k == 1) {
    int ibin = (nbinz-1)*nbinx*nbiny + j*nbinx + i;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,i,j,nbinz-1);
  }
  if (k == nbinz-2) {
    int ibin = j*nbinx + i;
    for (n = 0; n < nbinimages[m]; n++)
      if (binimage[m][n] == ibin) break;
    if (n == nbinimages[m]) binimage[m][nbinimages[m]++] = ibin;
    add_image_bins(m,i,j,0);
  }
}

/* ----------------------------------------------------------------------
   grow per-site arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AppOffLattice::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;

  id = (int *)  memory->srealloc(id,nmax*sizeof(int),"app:id");
  bin = (int *) memory->srealloc(bin,nmax*sizeof(int),"app:bin");
  next = (int *) memory->srealloc(next,nmax*sizeof(int),"app:next");
  prev = (int *) memory->srealloc(prev,nmax*sizeof(int),"app:prev");
  nextimage = (int *) memory->srealloc(nextimage,nmax*sizeof(int),
				       "app:nextimage");

  xyz = memory->grow_2d_double_array(xyz,nmax,3,"app:xyz");

  if (sitecustom == 0)
    site = (int *) memory->srealloc(site,nmax*sizeof(int),"app:site");
  else {
    for (int i = 0; i < ninteger; i++)
      iarray[i] = (int *) 
	memory->srealloc(iarray[i],nmax*sizeof(int),"app:iarray");
    for (int i = 0; i < ndouble; i++)
      darray[i] = (double *) 
	memory->srealloc(darray[i],nmax*sizeof(double),"app:darray");
  }
}

