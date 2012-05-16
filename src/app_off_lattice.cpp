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
#include "string.h"
#include "stdlib.h"
#include "app_off_lattice.h"
#include "comm_off_lattice.h"
#include "solve.h"
#include "domain.h"
#include "random_mars.h"
#include "random_park.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define DELTA 10000
#define MAXNEIGH 1000
#define EPSILON 1.0e-8

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
  site2i = NULL;
  in_sector = NULL;

  temperature = 0.0;

  nlocal = nghost = nmax = 0;
  nfree = 0;
  freehead = -1;

  bin = NULL;
  next = prev = nextimage = NULL;

  comm = NULL;

  binhead = binflag = NULL;
  pbcoffset = NULL;
  nimages = NULL;
  imageindex = imageproc = NULL;
  ghostindex = ghostproc = NULL;

  nstencil = 0;
  stencil = NULL;

  neighs = new int[MAXNEIGH];

  naccept = nattempt = 0;
  nsweeps = 0;
}

/* ---------------------------------------------------------------------- */

AppOffLattice::~AppOffLattice()
{
  delete [] set;

  delete ranapp;
  memory->destroy(sitelist);
  memory->destroy(site2i);
  memory->destroy(in_sector);

  memory->destroy(bin);
  memory->destroy(next);
  memory->destroy(prev);
  memory->destroy(nextimage);

  delete comm;

  memory->destroy(binhead);
  memory->destroy(binflag);
  memory->destroy(pbcoffset);
  memory->destroy(nimages);
  memory->destroy(imageindex);
  memory->destroy(imageproc);
  memory->destroy(ghostindex);
  memory->destroy(ghostproc);
  
  delete [] stencil;
  delete [] neighs;
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"sweep") == 0) set_sweep(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::init()
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

  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all(FLERR,"Cannot use KMC solver in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
    error->all(FLERR,
	       "Cannot use random rejection KMC in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
    error->all(FLERR,
	       "Cannot use raster rejection KMC in parallel with no sectors");

  if (sweepflag && dt_sweep == 0.0) 
    error->all(FLERR,"App did not set dt_sweep");

  // total size of a site = id + xyz + per-site quantities

  size_one = 4 + ninteger + ndouble;

  // domain settings

  dimension = domain->dimension;
  xprd = domain->xprd;
  yprd = domain->yprd;
  zprd = domain->zprd;
  subxlo = domain->subxlo;
  subylo = domain->subylo;
  subzlo = domain->subzlo;
  subxhi = domain->subxhi;
  subyhi = domain->subyhi;
  subzhi = domain->subzhi;

  // if sectors, set number of sectors

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

  // create sets based on sectors

  delete [] set;

  if (nsector == 1) {
    nset = 1;
    set = new Set[nset];
    create_set(0,0);
  } else {
    nset = nsector;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) create_set(i,i+1);
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
  // redo this on every init in case cutoff changed
  
  init_bins();

  // create stencil for neighbor finding
  // redo on every init in case cutoff and bins changed

  init_stencil();

  // bin owned sites
  // loop in reverse order so linked list will be in forward order

  for (int i = nlocal-1; i >= 0; i--) {
    bin[i] = site2bin(i);
    add_to_bin(i,bin[i]);
  }

  // initialize comm, both for this proc's full domain and sectors
  // redo on every init in case cutoff and bins changed

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

  // NOTE: dt_rkmc needs to be set appropriately
  // NOTE: test for bad choice of sector stop needs to be done

  dt_rkmc = dt_sweep;

  /*
  if (sweepflag != NOSWEEP) {
    double nme = 0.0;
    for (int i = 0; i < nset; i++) nme += set[i].nselect;
    double ntotal;
    MPI_Allreduce(&nme,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);

    dt_rkmc = ntotal/nglobal * dt_sweep;
    if (dt_rkmc == 0.0)
      error->all(FLERR,"Choice of sector stop led to no rKMC events");
    dt_rkmc = MIN(dt_rkmc,stoptime-time);
  }
  */

  // NOTE: sitelist needs to be dynamically allocated inside iter loop
  // ditto for site2i,etc ??  since nlocal can change dynamically

  // setup sitelist if sweepflag = RANDOM
  // do this every run since sector timestep could have changed

  if (sweepflag == RANDOM) {
    memory->destroy(sitelist);
    memory->create(sitelist,nlocal,"app:sitelist");
  }

  // setup site2i and in_sector lists
  // do this every run in case nlocal has changed

  memory->destroy(site2i);
  memory->destroy(in_sector);
  memory->create(site2i,nlocal,"app:site2i");
  memory->create(in_sector,nlocal,"app:in_sector");

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
  /*
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
  */
}

/* ----------------------------------------------------------------------
   KMC solver on sectors
   can be invoked in serial or parallel
 ------------------------------------------------------------------------- */

void AppOffLattice::iterate_kmc_sector(double stoptime)
{
  /*
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
  */
}

/* ----------------------------------------------------------------------
   rejection KMC solver
 ------------------------------------------------------------------------- */

void AppOffLattice::iterate_rejection(double stoptime)
{
  int i,nloop,nselect,nrange;

  // setup in_sector and site2i as if 1 processor with no sectors
  // will be reset in loop if there are sectors

  int nlocal_sector = nlocal;
  for (i = 0; i < nlocal; i++) {
    in_sector[i] = 1;
    site2i[i] = i;
  }
  
  int done = 0;
  while (!done) {
    for (int iset = 0; iset < nset; iset++) {
      activeset = iset;

      // acquire ghost sites for current sector

      //check("Before forward",0,iset);

      if (sectorflag) {
	timer->stamp();
	comm->sector(iset);
	timer->stamp(TIME_COMM);
      }

      //check("After forward",1,iset);

      timer->stamp();

      // determine which sites are currently in sector

      if (sectorflag) {
	nlocal_sector = 0;
	for (i = 0; i < nlocal; i++) {
	  in_sector[i] = inside_sector(i);
	  if (in_sector[i]) site2i[nlocal_sector++] = i;
	}
      }

      rkmc_params(nlocal_sector,nloop,nselect);

      // sweep over sites in sector (could be no sectors)
      // skip site if not currently in sector

      // random selection of sites in iset

      if (sweepflag == RANDOM) {
	nrange = nlocal_sector;
	for (i = 0; i < nselect; i++) 
	  sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
	for (int m = 0; m < nselect; m++) {
	  if (in_sector[sitelist[m]] == 0) continue;
	  site_event_rejection(sitelist[m],ranapp);
	}
	nattempt += nselect;

      // ordered sweep over all sites in iset

      } else {
	for (i = 0; i < nloop; i++)
	  for (int m = 0; m < nlocal_sector; m++) {
	    if (i && in_sector[site2i[m]] == 0) continue;
	    site_event_rejection(site2i[m],ranapp);
	  }
	nattempt += nselect;
      }

      timer->stamp(TIME_SOLVE);

      // migrate sites that moved into ghost bins

      //check("Before reverse",1,iset);

      if (sectorflag) {
	comm->reverse_sector(iset);
	timer->stamp(TIME_COMM);
      }

      //check("After reverse",0,iset);
    }

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ----------------------------------------------------------------------
  convert rejection info to rKMC params
 ------------------------------------------------------------------------- */

void AppOffLattice::rkmc_params(int nlocal_sector, int &nloop, int &nselect)
{
  if (sweepflag == RANDOM) {
    if (nstop > 0.0) {
      nloop = 0;
      nselect = static_cast<int> (nstop*nlocal_sector);
    }
    if (tstop > 0.0) {
      double n = tstop / (dt_sweep/nglobal);
      nloop = 0;
      nselect = static_cast<int> (n/nglobal * nlocal_sector);
    }
    
  } else if (sweepflag == RASTER) {
    int n;
    if (nstop > 0.0) n = static_cast<int> (nstop);
    if (tstop > 0.0) n = static_cast<int> (tstop/dt_sweep);
    nloop = n;
    nselect = n * nlocal_sector;
  }
}

/* ----------------------------------------------------------------------
   debug check that bin data structures are self-consistent
 ------------------------------------------------------------------------- */

void AppOffLattice::check(char *str, int flag, int iset)
{
  int nall = nlocal;
  if (flag) nall = nlocal+nghost;

  // set pflag to 1 if want success output every step
  // set maxcount to reasonable number of sites per bin

  int pflag = 0;
  int maxcount = 40;

  MPI_Barrier(world);

  // test that bin index of every site is correct and within bin bounds

  for (int i = 0; i < nall; i++) {
    if (site2bin(i) != bin[i] || bin[i] < 0 || bin[i] > nbins-1) {
      printf("%s SITE %d: %d " TAGINT_FORMAT " %d: %g %g %g\n",
	     str,me,i,id[i],bin[i],xyz[i][0],xyz[i][1],xyz[i][2]);
      error->one(FLERR,"SITE MISMATCH");
    }
  }

  // test that owned bins only contain owned sites

  for (int i = 0; i < nbins; i++) {
    if (binflag[i] == GHOST) continue;
    int m = binhead[i];
    while (m >= 0) {
      if (m >= nlocal) {
	printf("%s GHOST %d: %d %d %d\n",str,me,m,nlocal,i);
	error->one(FLERR,"GHOST IN OWNED BIN");
      }
      m = next[m];
    }
  }

  // test that each bin's linked list is consistent

  for (int i = 0; i < nbins; i++) {
    int m = binhead[i];
    int mprev = -1;
    while (m >= 0) {
      if (prev[m] != mprev) {
	printf("%s LINK me %d: m %d id " TAGINT_FORMAT 
	       " bin %d prev %d mprev %d nloc %d\n",
	       str,me,m,id[m],bin[m],prev[m],mprev,nlocal);
	error->one(FLERR,"LINK MISMATCH");
      }
      mprev = m;
      m = next[m];
    }
  }

  // test that site bin matches bin linked list it is in

  for (int i = 0; i < nbins; i++) {
    int m = binhead[i];
    while (m >= 0) {
      if (bin[m] != i) {
	printf("%s BIN %d: %d " TAGINT_FORMAT " %d %d %d\n",
	       str,me,m,id[m],bin[m],i,site2bin(m));
	error->one(FLERR,"BIN MISMATCH");
      }
      m = next[m];
    }
  }

  // test that walking lists of all bins sums to total site count

  int count = 0;
  for (int i = 0; i < nbins; i++) {
    if (!flag && binflag[i] == GHOST) continue;
    int m = binhead[i];
    while (m >= 0) {
      count++;
      m = next[m];
    }
  }
  if (count != nall) {
    printf("SITES NOT IN BINS %d %d %d %d %d\n",
	   count,nlocal,nghost,nall,flag);
    error->one(FLERR,"SITES NOT IN BINS");
  }

  // test that no bin has too many sites
  
  for (int i = 0; i < nbins; i++) {
    count = 0;
    int m = binhead[i];
    while (m >= 0) {
      count++;
      m = next[m];
    }
    if (count > maxcount) {
      printf("%s COUNT %d %d %d\n",str,me,i,count);
      error->one(FLERR,"COUNT MISMATCH");
    }
  }

  if (pflag && me == 0) printf("OK %d %g %s %d\n",me,time,str,iset);
  MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::input_app(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::set_sector(int narg, char **arg)
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

void AppOffLattice::set_sweep(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal sweep command");
  if (strcmp(arg[0],"random") == 0) sweepflag = RANDOM;
  else if (strcmp(arg[0],"raster") == 0) sweepflag = RASTER;
  else if (strcmp(arg[0],"none") == 0) sweepflag = NOSWEEP;
  else error->all(FLERR,"Illegal sweep command");
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppOffLattice::stats(char *strtmp)
{
  char big[8],format[64];
  strcpy(big,BIGINT_FORMAT);

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  if (solve) {
    sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
    sprintf(strtmp,format,time,naccept_all,0,0);
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

void AppOffLattice::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nreject","Nsweeps");
}

/* ----------------------------------------------------------------------
   create a subset geometry within my sub-domain
   isector = 0 = all sites (no sector)
   isector >= 1 = sites within a sector
 ------------------------------------------------------------------------- */

void AppOffLattice::create_set(int iset, int isector)
{
  if (isector == 0) {
    set[iset].xlo = subxlo;
    set[iset].xhi = subxhi;
    set[iset].ylo = subylo;
    set[iset].yhi = subyhi;
    set[iset].zlo = subzlo;
    set[iset].zhi = subzhi;

  } else {
    int ix = (isector-1) % 2;
    int iy = ((isector-1)/2) % 2;
    int iz = (isector-1) / 4;

    if (ix == 0) {
      set[iset].xlo = subxlo;
      set[iset].xhi = 0.5 * (subxlo + subxhi);
    } else {
      set[iset].xlo = 0.5 * (subxlo + subxhi);
      set[iset].xhi = subxhi;
    }

    if (dimension == 1) {
      set[iset].ylo = subylo;
      set[iset].yhi = subyhi;
    } else if (iy == 0) {
      set[iset].ylo = subylo;
      set[iset].yhi = 0.5 * (subylo + subyhi);
    } else {
      set[iset].ylo = 0.5 * (subylo + subyhi);
      set[iset].yhi = subyhi;
    }

    if (dimension <= 2) {
      set[iset].zlo = subzlo;
      set[iset].zhi = subzhi;
    } else if (iz == 0) {
      set[iset].zlo = subzlo;
      set[iset].zhi = 0.5 * (subzlo + subzhi);
    } else {
      set[iset].zlo = 0.5 * (subzlo + subzhi);
      set[iset].zhi = subzhi;
    }
  }
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::init_bins()
{
  int i,j,k,m,ix,iy,iz;

  int *procgrid = domain->procgrid;
  int *myloc = domain->myloc;

  // bin size must be >= cutoff set by parameters in init_app()
  // require 2 bins within sub-domain if sectoring in a dimension
  // else require 1 bin within sub-domain
  // require an even number of bins if sectoring in a dimension

  double cutoff = delpropensity + delevent;
  nbinx = static_cast<int> (xprd/procgrid[0]/cutoff);
  nbiny = nbinz = 1;
  if (dimension >= 2) nbiny = static_cast<int> (yprd/procgrid[1]/cutoff);
  if (dimension == 3) nbinz = static_cast<int> (zprd/procgrid[2]/cutoff);

  if (nsector >= 2) nbinx = nbinx/2 * 2;
  if (nsector >= 4) nbiny = nbiny/2 * 2;
  if (nsector == 8) nbinz = nbinz/2 * 2;

  if (nsector >= 2 && nbinx < 2)
    error->all(FLERR,"Application cutoff is too big for processor sub-domain");
  if (nsector >= 4 && nbiny < 2)
    error->all(FLERR,"Application cutoff is too big for processor sub-domain");
  if (nsector == 8 && nbinz < 2)
    error->all(FLERR,"Application cutoff is too big for processor sub-domain");
  if (nbinx < 1 || nbiny < 1 || nbinz < 1)
    error->all(FLERR,"Application cutoff is too big for processor sub-domain");

  binx = xprd/procgrid[0] / nbinx;
  biny = yprd/procgrid[1] / nbiny;
  binz = zprd/procgrid[2] / nbinz;
  invbinx = 1.0 / binx;
  invbiny = 1.0 / biny;
  invbinz = 1.0 / binz;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Bins/proc: %d %d %d\n",nbinx,nbiny,nbinz);
      fprintf(screen,"Bin sizes: %g %g %g\n",binx,biny,binz);
    }
    if (logfile) {
      fprintf(logfile,"Bins/proc: %d %d %d\n",nbinx,nbiny,nbinz);
      fprintf(logfile,"Bin sizes: %g %g %g\n",binx,biny,binz);
    }
  }

  // add 2 bins for ghost regions

  nbinx += 2;
  if (dimension >= 2) nbiny += 2;
  if (dimension == 3) nbinz += 2;

  nbins = nbinx*nbiny*nbinz;

  memory->destroy(binhead);
  memory->destroy(binflag);
  memory->destroy(pbcoffset);
  memory->destroy(nimages);
  memory->destroy(imageindex);
  memory->destroy(ghostindex);
  memory->destroy(ghostproc);

  memory->create(binhead,nbins,"app:binhead");
  memory->create(binflag,nbins,"app:binflag");
  memory->create(pbcoffset,nbins,3,"app:pbcoffset");
  memory->create(nimages,nbins,"app:nimages");
  if (dimension == 3) {
    memory->create(imageindex,nbins,7,"app:imageindex");
    memory->create(imageproc,nbins,7,"app:imageproc");
  } else if (dimension == 2) {
    memory->create(imageindex,nbins,3,"app:imageindex");
    memory->create(imageproc,nbins,3,"app:imageproc");
  } else {
    memory->create(imageindex,nbins,1,"app:imageindex");
    memory->create(imageproc,nbins,1,"app:imageproc");
  }
  memory->create(ghostindex,nbins,"app:ghostindex");
  memory->create(ghostproc,nbins,"app:ghostproc");

  // determine if each bin is INTERIOR or EDGE or GHOST

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

  // set binhead to empty for all bins
  // set pbcoffset, ghostindex, ghostproc, for each GHOST bin
  // set nimages, imageindex, imageproc for each EDGE bin

  for (m = 0; m < nbins; m++) {
    binhead[m] = -1;
    if (binflag[m] == GHOST) {
      ix = m % nbinx;
      iy = (m/nbinx) % nbiny;
      iz = m / (nbinx*nbiny);

      if (ix == 0 && myloc[0] == 0) pbcoffset[m][0] = -1;
      else if (ix == nbinx-1 && myloc[0] == procgrid[0]-1)
	pbcoffset[m][0] = 1;
      else pbcoffset[m][0] = 0;

      if (dimension >= 2) {
	if (iy == 0 && myloc[1] == 0) pbcoffset[m][1] = -1;
	else if (iy == nbiny-1 && myloc[1] == procgrid[1]-1) 
	  pbcoffset[m][1] = 1;
	else pbcoffset[m][1] = 0;
      } else pbcoffset[m][1] = 0;

      if (dimension == 3) {
	if (iz == 0 && myloc[2] == 0) pbcoffset[m][2] = -1;
	else if (iz == nbinz-1 && myloc[2] == procgrid[2]-1) 
	  pbcoffset[m][2] = 1;
	else pbcoffset[m][2] = 0;
      } else pbcoffset[m][2] = 0;

      ghostproc[m] = neighproc(1,ix,iy,iz);

      if (ix == 0) ix = nbinx-2;
      if (ix == nbinx-1) ix = 1;

      if (dimension >= 2) {
	if (iy == 0) iy = nbiny-2;
	if (iy == nbiny-1) iy = 1;
      } else iy = 0;

      if (dimension == 3) {
	if (iz == 0) iz = nbinz-2;
	if (iz == nbinz-1) iz = 1;
      } else iz = 0;

      ghostindex[m] = iz*nbinx*nbiny + iy*nbinx + ix;

    } else if (binflag[m] == EDGE) {
      ix = m % nbinx;
      iy = (m/nbinx) % nbiny;
      iz = m / (nbinx*nbiny);
      nimages[m] = 0;
      add_image_bins(m,ix,iy,iz);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AppOffLattice::init_stencil()
{
  int i,j,k;

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

    // DEBUG CHECK - can remove
    if (ibin+stencil[k] < 0 || ibin+stencil[k] >= nbins) {
      printf("BAD STENCIL %d %d %d\n",i,ibin,ibin+stencil[k]);
      error->one(FLERR,"BAD STENCIL");
    }

    for (j = binhead[ibin+stencil[k]]; j >= 0; j = next[j]) {
      delx = xtmp - xyz[j][0];
      dely = ytmp - xyz[j][1];
      delz = ztmp - xyz[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // SHOULD DELETE THIS TEST eventually
      if (numneigh == MAXNEIGH) {
	printf("NEIGH %d %d %d\n",i,ibin,numneigh);
	error->one(FLERR,"Too many neighbors per site");
      }

      if (rsq <= cutsq && i != j) neighs[numneigh++] = j;
    }
  }

  if (numneigh > MAXNEIGH) error->one(FLERR,"Too many neighbors per site");
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

  // if sectoring (1 or more procs), then move is simple
  // if stays in same bin, done
  // else delete it from old bin, add it to new bin
  // also flag site if moves out of sector

  if (sectorflag) {
    if (newbin == oldbin) return;
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    in_sector[i] = inside_sector(i);
    return;
  }

  // remaining logic is all for 1 processor with no sectors
  // complicated due to ghost images and PBC
  // site is remapped for PBC, so is never flagged as out of sector
  // 6 different cases to consider: A-F

  // (A) site stays in same bin
  // update its ghost image coords if necessary

  if (newbin == oldbin) {
    j = nextimage[i];
    while (j >= 0) {
      oldbin = bin[j];
      xyz[j][0] = xyz[i][0] + pbcoffset[oldbin][0]*xprd;
      xyz[j][1] = xyz[i][1] + pbcoffset[oldbin][1]*yprd;
      xyz[j][2] = xyz[i][2] + pbcoffset[oldbin][2]*zprd;
      j = nextimage[j];
    }
    return;
  }

  // (B) site stays in INTERIOR
  // delete it from old bin, add it to new bin

  if (binflag[oldbin] == INTERIOR && binflag[newbin] == INTERIOR) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    return;
  }

  // (C) site stays in EDGE
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

  // (D) site moves from INTERIOR to EDGE
  // delete it from old bin, add it to new bin
  // create new ghost images

  if (binflag[oldbin] == INTERIOR && binflag[newbin] == EDGE) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    add_images(i,newbin);
    return;
  }

  // (E) site moves from EDGE to INTERIOR
  // delete it from old bin, add it to new bin
  // delete ghost images

  if (binflag[oldbin] == EDGE && binflag[newbin] == INTERIOR) {
    delete_from_bin(i,oldbin);
    add_to_bin(i,newbin);
    bin[i] = newbin;
    delete_images(i);
    return;
  }

  // (F) site moves from EDGE to GHOST
  // map it back into simulation box via PBC, then is same as EDGE to EDGE
  // delete it from old bin, add it to new bin
  // delete ghost images
  // create new ghost images

  if (binflag[oldbin] == EDGE && binflag[newbin] == GHOST) {
    xyz[i][0] -= pbcoffset[newbin][0]*xprd;
    xyz[i][1] -= pbcoffset[newbin][1]*yprd;
    xyz[i][2] -= pbcoffset[newbin][2]*zprd;
    newbin = site2bin(i);
    if (binflag[newbin] != EDGE) error->one(FLERR,"PBC remap of site failed");
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
  // DEBUG check - can delete eventually

  if (xyz[i][0] < subxlo-binx-EPSILON || 
      xyz[i][0] >= subxhi+binx+EPSILON ||
      xyz[i][1] < subylo-biny-EPSILON || 
      xyz[i][1] >= subyhi+biny+EPSILON ||
      xyz[i][2] < subzlo-binz-EPSILON || 
      xyz[i][2] >= subzhi+binz+EPSILON) {
    printf("BAD SITE: %d %d %d %d: %g %g %g " TAGINT_FORMAT "\n",
	   me,i,bin[i],nlocal,xyz[i][0],xyz[i][1],xyz[i][2],id[i]);
    printf("MY DOMAIN: %g %g %g: %g %g %g\n",subxlo,subylo,
	   subzlo,subxhi,subyhi,subzhi);
    error->one(FLERR,"Site not in my bin domain");
  }

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
   add site I to beginning of bin IBIN
 ------------------------------------------------------------------------- */

void AppOffLattice::add_to_bin(int i, int ibin)
{
  // DEBUG check that I is in IBIN - can delete eventually

  if (ibin < 0 || ibin >= nbins)
    error->one(FLERR,"Adding site to illegal bin");

  if (site2bin(i) != ibin) {
    printf("BAD BIN: %d %d %d: %g %g %g\n",
	   me,i,nlocal,xyz[i][0],xyz[i][1],xyz[i][2]);
    printf("IBIN: %d\n",ibin);
    error->one(FLERR,"Adding site to bin it is not in");
  }

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
   only called when 1 proc and no sectors, so all images are owned by proc 0
 ------------------------------------------------------------------------- */

void AppOffLattice::add_images(int i, int ibin)
{
  int j,jbin;

  nextimage[i] = -1;
  int n = nimages[ibin];
  for (int m = 0; m < n; m++) {
    jbin = imageindex[ibin][m];
    j = new_ghost_site();
    xyz[j][0] = xyz[i][0] + pbcoffset[jbin][0]*xprd;
    xyz[j][1] = xyz[i][1] + pbcoffset[jbin][1]*yprd;
    xyz[j][2] = xyz[i][2] + pbcoffset[jbin][2]*zprd;
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

int AppOffLattice::new_ghost_site()
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
   return index of a new owned site from free list
   grow per-site arrays if free list is empty
 ------------------------------------------------------------------------- */

int AppOffLattice::new_owned_site()
{
  if (nfree == 0) {
    int oldmax = nmax;
    grow(0);
    add_free(oldmax);
  }

  int index = freehead;
  freehead = next[freehead];
  nfree--;
  nlocal++;
  return index;
}

/* ----------------------------------------------------------------------
   delete an owned site
   keep owned list compact by copying last site into deleted site
   add last site to free list
   return next ptr of deleted site
 ------------------------------------------------------------------------- */

int AppOffLattice::delete_owned_site(int i)
{
  int nextptr = next[i];
  if (nextptr == nlocal-1) nextptr = i;

  delete_from_bin(i,bin[i]);

  id[i] = id[nlocal-1];
  xyz[i][0] = xyz[nlocal-1][0];
  xyz[i][1] = xyz[nlocal-1][1];
  xyz[i][2] = xyz[nlocal-1][2];
  for (int k = 0; k < ninteger; k++) iarray[k][i] = iarray[k][nlocal-1];
  for (int k = 0; k < ndouble; k++) darray[k][i] = darray[k][nlocal-1];

  bin[i] = bin[nlocal-1];
  next[i] = next[nlocal-1];
  prev[i] = prev[nlocal-1];
  nextimage[i] = nextimage[nlocal-1];

  // reset next/prev ptrs of adjacent linked sites
  // do not call delete_from_bin() and add_to_bin() because
  //   are looping over sites in I's bin,
  //   and loop will be messed up if nlocal-1 site is also in the bin

  if (prev[i] < 0) binhead[bin[i]] = i;
  else next[prev[i]] = i;
  if (next[i] >= 0) prev[next[i]] = i;

  next[nlocal-1] = freehead;
  freehead = nlocal-1;
  nfree++;
  nlocal--;

  return nextptr;
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
   generate image bins of edge bin M, which is at local indices I,J,K
   called recursively with ghost bins to generate all images
   image = same bin as M, stored as ghost bin by some proc
   imageproc = proc that owns ghost bin, could be self
   imageindex = local index of ghost bin on owning proc
   inew,jnew,knew = ghost bin indices on proc that owns it
   do not add image bin if already in image[] list
   can add up to 3 image bins in 2d, 7 in 3d
------------------------------------------------------------------------- */

void AppOffLattice::add_image_bins(int m, int i, int j, int k)
{
  int n,inew,jnew,knew,newbin;

  if (i == 1) {
    inew = nbinx-1;
    jnew = j;
    knew = k;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
  }
  if (i == nbinx-2) {
    inew = 0;
    jnew = j;
    knew = k;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
  }

  if (j == 1) {
    inew = i;
    jnew = nbiny-1;
    knew = k;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
  }
  if (j == nbiny-2) {
    inew = i;
    jnew = 0;
    knew = k;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
  }

  if (k == 1) {
    inew = i;
    jnew = j;
    knew = nbinz-1;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
  }
  if (k == nbinz-2) {
    inew = i;
    jnew = j;
    knew = 0;
    newbin = knew*nbinx*nbiny + jnew*nbinx + inew;
    for (n = 0; n < nimages[m]; n++)
      if (imageindex[m][n] == newbin) break;
    if (n == nimages[m]) {
      imageindex[m][nimages[m]] = newbin;
      imageproc[m][nimages[m]] = neighproc(-1,inew,jnew,knew);
      nimages[m]++;
    }
    add_image_bins(m,inew,jnew,knew);
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
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  memory->grow(id,nmax,"app:id");
  memory->grow(bin,nmax,"app:bin");
  memory->grow(next,nmax,"app:next");
  memory->grow(prev,nmax,"app:prev");
  memory->grow(nextimage,nmax,"app:nextimage");

  memory->grow(xyz,nmax,3,"app:xyz");

  memory->grow(site2i,nmax,"app:site2i");
  memory->grow(in_sector,nmax,"app:in_sector");

  for (int i = 0; i < ninteger; i++)
    memory->grow(iarray[i],nmax,"app:iarray");
  for (int i = 0; i < ndouble; i++)
    memory->srealloc(darray[i],nmax,"app:darray");

  grow_app();
}

/* ----------------------------------------------------------------------
   return proc ID of owner of ghost bin with local indices inew,jnew,knew
   returned proc ID could be self due to PBC
------------------------------------------------------------------------- */

int AppOffLattice::neighproc(int offset, int inew, int jnew, int knew)
{
  int i,j,k,idelta,jdelta,kdelta;

  // i,j,k = indices of me in 3d grid of procs

  int *procgrid = domain->procgrid;

  i = me % procgrid[0];
  j = (me/procgrid[0]) % procgrid[1];
  k = me / (procgrid[0]*procgrid[1]);

  // idelta,jdelta,kdelta = which of my neighbors the ghost bin maps to

  if (inew == nbinx-1) idelta = offset;
  else if (inew == 0) idelta = -offset;
  else idelta = 0;
  if (jnew == nbiny-1) jdelta = offset;
  else if (jnew == 0) jdelta = -offset;
  else jdelta = 0;
  if (knew == nbinz-1) kdelta = offset;
  else if (knew == 0) kdelta = -offset;
  else kdelta = 0;

  if (dimension == 1) jdelta = kdelta = 0;
  if (dimension == 2) kdelta = 0;

  // i,j,k = indices of neighbor proc in 3d grid of procs

  i += idelta;
  if (i < 0) i = procgrid[0]-1;
  if (i == procgrid[0]) i = 0;
  j += jdelta;
  if (j < 0) j = procgrid[1]-1;
  if (j == procgrid[1]) j = 0;
  k += kdelta;
  if (k < 0) k = procgrid[2]-1;
  if (k == procgrid[2]) k = 0;

  int newproc = k*procgrid[0]*procgrid[1] + j*procgrid[0] + i;
  return newproc;
}

/* ----------------------------------------------------------------------
   return 1 if site I is inside activeset sector, else return 0
   activeset could be entire proc sub-domain (no sectors)
------------------------------------------------------------------------- */

int AppOffLattice::inside_sector(int i)
{
  if (xyz[i][0] < set[activeset].xlo) return 0;
  if (xyz[i][0] >= set[activeset].xhi) return 0;
  if (xyz[i][1] < set[activeset].ylo) return 0;
  if (xyz[i][1] >= set[activeset].yhi) return 0;
  if (xyz[i][2] < set[activeset].zlo) return 0;
  if (xyz[i][2] >= set[activeset].zhi) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   add an owned site
   called from create_sites or read_sites commands
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void AppOffLattice::add_site(tagint n, double x, double y, double z)
{
  if (nlocal == nmax) grow(0);

  id[nlocal] = n;
  xyz[nlocal][0] = x;
  xyz[nlocal][1] = y;
  xyz[nlocal][2] = z;

  for (int i = 0; i < ninteger; i++) iarray[i][nlocal] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][nlocal] = 0;

  nlocal++;
}

/* ----------------------------------------------------------------------
   set values for owned site I
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppOffLattice::add_values(int i, char **values)
{
  for (int m = 0; m < ninteger; m++) iarray[m][i] = atoi(values[m]);
  for (int m = 0; m < ndouble; m++) darray[m][i] = atof(values[m+ninteger]);
}
