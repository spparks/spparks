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
#include "app_lattice2d.h"
#include "sweep_lattice2d.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "output.h"

using namespace SPPARKS_NS;

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

AppLattice2d::AppLattice2d(SPPARKS *spk, int narg, char **arg) : 
  App(spk,narg,arg)
{
  appclass = LATTICE2D;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  ntimestep = 0;
  time = 0.0;
  temperature = 0.0;
  ibufread = NULL;
  fp = NULL;
  propensity = NULL;
  site2ij = NULL;
  ij2site = NULL;
  id = NULL;
  xyz = NULL;

  // setup communicator for ghost sites

  comm = new CommLattice2d(spk);

  // app can override these values in its constructor

  delpropensity = 1;
  delevent = 0;
  numrandom = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice2d::~AppLattice2d()
{
  delete comm;

  delete [] ibufread;
  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ij);
  memory->destroy_2d_T_array(ij2site);
  memory->sfree(id);
  memory->destroy_2d_T_array(xyz);
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::init()
{
  int i,j,m;

  // error checks

  if (sweep && strcmp(sweep->style,"lattice2d") != 0)
    error->all("Mismatched sweeper with app lattice");

  if (sweep == NULL && solve == NULL)
    error->all("Lattice app needs a solver or sweeper");

  if (sweep && ((SweepLattice2d *) sweep)->Lkmc && solve == NULL)
    error->all("Must define solver with KMC sweeper");

  if (solve && sweep && ((SweepLattice2d *) sweep)->Lkmc == false)
    error->all("Cannot use solver with non-KMC sweeper");

  if (solve && sweep == NULL && nprocs > 1)
    error->all("Cannot use solver in parallel");

  // app-specific initialization

  init_app();

  // comm init

  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth,
	     delpropensity,delevent);

  // if no sweeper, initialize 3 arrays: propensity, site2i,i2site
  // sweeper allocates its own per-sector versions of these

  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ij);
  memory->destroy_2d_T_array(ij2site);

  int nsites = nx_local*ny_local;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "app2d:propensity");
    memory->create_2d_T_array(site2ij,nsites,2,
			      "app2d:site2ij");
    memory->create_2d_T_array(ij2site,nxlo,nxhi,nylo,nyhi,
			      "app2d:ij2site");

    for (m = 0; m < nsites; m++) {
      i = m / ny_local + 1;
      j = m % ny_local + 1;
      site2ij[m][0] = i;
      site2ij[m][1] = j;
    }

    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++)
	ij2site[i][j] = (i-1)*ny_local + j-1;

  } else {
    propensity = NULL;
    site2ij = NULL;
    ij2site = NULL;
  }

  // initialize sweeper
  
  if (sweep) {
    sweep->init();
    Lmask = sweep->Lmask;
    mask = ((SweepLattice2d *) sweep)->mask;
  }

  // initialize propensities for KMC solver
  // if KMC sweep, sweeper does its own init of its propensity arrays
  // comm insures ghost sites are set

  if (sweep == NULL) {
    comm->all(lattice);
    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++)
	propensity[ij2site[i][j]] = site_propensity(i,j);
    solve->init(nsites,propensity);
  }

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) output->set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::input_app(char *command, int narg, char **arg)
{
  error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice2d::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

  init();
  timer->init();

  iterate();
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on solver
 ------------------------------------------------------------------------- */

void AppLattice2d::iterate()
{
  int i,j,isite;
  double dt;
    
  if (time >= stoptime) return;
  timer->barrier_start(TIME_LOOP);
  
  int done = 0;
  while (!done) {
    if (sweep) {
      sweep->do_sweep(dt);
      time += dt;

    } else {
      timer->stamp();
      isite = solve->event(&dt);
      timer->stamp(TIME_SOLVE);
      
      if (isite < 0) done = 1;
      else {
	ntimestep++;
	i = site2ij[isite][0];
	j = site2ij[isite][1];
	site_event(i,j,1,random);
	time += dt;
	timer->stamp(TIME_APP);
      }
    }

    if (time >= stoptime) done = 1;

    output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   update all ghost images of site i,j
   called when performing serial KMC on entire domain
------------------------------------------------------------------------- */

void AppLattice2d::update_ghost_sites(int i, int j)
{
  // i = 1 line becomes i = nx_local+1 line w/ j ghosts
  // i = nx_local line becomes i = 0 line w/ j ghosts

  if (i == 1) {
    lattice[nx_local+1][j] = lattice[i][j];
    if (j == 1) lattice[nx_local+1][ny_local+1] = lattice[i][j];
    if (j == ny_local) lattice[nx_local+1][0] = lattice[i][j];
  }
  if (i == nx_local) {
    lattice[0][j] = lattice[i][j];
    if (j == 1) lattice[0][ny_local+1] = lattice[i][j];
    if (j == ny_local) lattice[0][0] = lattice[i][j];
  }

  // j = 1 line becomes j = ny_local+1 line w/out i ghosts
  // j = ny_local line becomes j = 0 line w/out i ghosts

  if (j == 1) lattice[i][ny_local+1] = lattice[i][j];
  if (j == ny_local) lattice[i][0] = lattice[i][j];
}

/* ----------------------------------------------------------------------
   add i,j value to site list if different than oldstate and existing sites
------------------------------------------------------------------------- */

void AppLattice2d::add_unique(int oldstate, int &nevent, int *sites,
			      int i, int j)
{
  int value = lattice[i][j];
  if (value == oldstate) return;
  for (int m = 0; m < nevent; m++)
    if (value == sites[m]) return;
  sites[nevent++] = value;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice2d::stats(char *strtmp)
{
  int ntimestepall;
  MPI_Allreduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,world);
  sprintf(strtmp," %10d %10g",ntimestepall,time);
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice2d::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s","Step","Time");
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d global lattice so as to minimize perimeter per proc
   setup id and xyz arrays
------------------------------------------------------------------------- */

void AppLattice2d::procs2lattice()
{
  int ipx,ipy;
  double boxx,boxy,surf;
  double bestsurf = 2 * (nx_global + ny_global);
  
  // loop thru all possible factorizations of nprocs
  // surf = perimeter of a proc sub-domain
 
  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = float(nx_global)/ipx;
      boxy = float(ny_global)/ipy;
      surf = boxx + boxy;
      if (surf < bestsurf) {
	bestsurf = surf;
	nx_procs = ipx;
	ny_procs = ipy;
      }
    }
    ipx++;
  }

  int iprocx = me/ny_procs;
  nx_offset = iprocx*nx_global/nx_procs;
  nx_local = (iprocx+1)*nx_global/nx_procs - nx_offset;

  int iprocy = me % ny_procs;
  ny_offset = iprocy*ny_global/ny_procs;
  ny_local = (iprocy+1)*ny_global/ny_procs - ny_offset;

  if (delevent > delpropensity) 
    error->all("Delevent > delpropensity");
  if (nx_local < 2*delpropensity || ny_local < 2*delpropensity)
    error->one("Lattice per proc is too small");

  if (iprocy == 0) procsouth = me + ny_procs - 1;
  else procsouth = me - 1;
  if (iprocy == ny_procs-1) procnorth = me - ny_procs + 1;
  else procnorth = me + 1;

  if (iprocx == 0) procwest = me + nprocs - ny_procs;
  else procwest = me - ny_procs;
  if (iprocx == nx_procs-1) proceast = me - nprocs + ny_procs;
  else proceast = me + ny_procs;

  nxlo = 1 - delpropensity;
  nxhi = nx_local + delpropensity;
  nylo = 1 - delpropensity;
  nyhi = ny_local + delpropensity;

  nglobal = nx_global * ny_global;
  nlocal = nx_local * ny_local;
  boxxlo = boxylo = 0.0;
  boxxhi = nx_global;
  boxyhi = ny_global;
  id = (int *) memory->smalloc(nlocal*sizeof(int),"app2d:id");
  memory->create_2d_T_array(xyz,nlocal,3,"app2d:xyz");

  int i,j;
  int m = 0;
  for (j = 0; j < ny_local; j++)
    for (i = 0; i < nx_local; i++) {
      id[m] = (j+ny_offset)*nx_global + i+nx_offset + 1;
      xyz[m][0] = i+nx_offset;
      xyz[m][1] = j+ny_offset;
      xyz[m][2] = 0.0;
      m++;
    }
}

/* ----------------------------------------------------------------------
   dump a snapshot of grains to the screen in 2D layout
   all the spins for each processor domain are printed out
   to call this function, follow this example:

//   char *fstring2 = "In stats() after comm->all() Timestep = %d";
//   int len2 = strlen(fstring2)+32;
//   char *title2 = new char[len2];
//   sprintf(title2,fstring2,ntimestep);
//   dump_detailed(title2);
//   delete [] title2;

------------------------------------------------------------------------- */

void AppLattice2d::dump_detailed(char* title)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxhtmp,nyhtmp,nzhtmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff

  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delpropensity)*
    ((ny_global-1)/ny_procs+1+2*delpropensity)+9;
  nsend = (nx_local+2*delpropensity)*(ny_local+2*delpropensity)+9;
  if (maxbuftmp < nsend) 
    error->one("Maxbuftmp size too small in AppGrain::dump_detailed()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),
				  "app2d:dump_detailed:buftmp");

  // proc 0 writes interactive dump header

  if (me == 0) {
    if (screen) {
      fprintf(screen,"*** Detailed Dump ***\n");
      fprintf(screen,"Title = %s\n",title);
      fprintf(screen,"nx_global = %d ny_global = %d\n",nx_global,ny_global);
    }
  }

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = 0;
  // Need to delete these two
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = 0;

  // pack my lattice values into buffer
  // Need to violate normal ordering in order to simplify output

  for (int j = 1-delpropensity; j <= ny_local+delpropensity; j++) {
    for (int i = 1-delpropensity; i <= nx_local+delpropensity; i++) {
      buftmp[m++] = lattice[i][j];
    }
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;

      if (screen) {
	m = 0;
	nxtmp = buftmp[m++];
	nytmp = buftmp[m++];
	nztmp = buftmp[m++];
	nxhtmp = buftmp[m++];
	nyhtmp = buftmp[m++];
	nzhtmp = buftmp[m++];
	nxotmp = buftmp[m++];
	nyotmp = buftmp[m++];
	nzotmp = buftmp[m++];
	fprintf(screen,"iproc = %d \n",iproc);
	fprintf(screen,"nxlocal = %d \nnylocal = %d \n",nxtmp,nytmp);
	fprintf(screen,"nx_offset = %d \nny_offset = %d \n",nxotmp,nyotmp);
	m = nrecv;
	for (int j = nytmp+delpropensity; j >= 1-delpropensity; j--) {
	  m-=nxtmp+2*delpropensity;
	  for (int i = 1-delpropensity; i <= nxtmp+delpropensity; i++) {
	    fprintf(screen,"%3d",buftmp[m++]);
	  }
	  fprintf(screen,"\n");
	  m-=nxtmp+2*delpropensity;
	}
      }
      
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
}

/* ----------------------------------------------------------------------
   dump a mask snapshot to the screen in 2D layout
   all the mask values for each processor domain are printed out
   to call this function, follow this example:

//       char *fstring1 = "Before sector() icolor = %d iquad = %d ";
//       int len1 = strlen(fstring1)+32;
//       char *title1 = new char[len1];
//       sprintf(title1,fstring1,icolor,iquad);
//       applattice->dump_detailed(title1);
//       applattice->dump_detailed_mask(title1,mask);
//       delete [] title1;

------------------------------------------------------------------------- */

void AppLattice2d::dump_detailed_mask(char* title, char** mask)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxhtmp,nyhtmp,nzhtmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff

  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delpropensity)*
    ((ny_global-1)/ny_procs+1+2*delpropensity)+9;
  nsend = (nx_local+2*delpropensity)*(ny_local+2*delpropensity)+9;
  if (maxbuftmp < nsend) 
    error->one("Maxbuftmp size too small in AppGrain::dump_detailed_mask()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),
				  "app2d:dump_detailed_mask:buftmp");

  // proc 0 writes interactive dump header

  if (me == 0) {
    if (screen) {
      fprintf(screen,"*** Detailed Mask Dump ***\n");
      fprintf(screen,"Title = %s\n",title);
      fprintf(screen,"nx_global = %d ny_global = %d\n",nx_global,ny_global);
    }
  }

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = 0;
  // Need to delete these two
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = 0;

  // pack the mask values into buffer
  // Need to violate normal ordering in order to simplify output

  for (int j = 1-delpropensity; j <= ny_local+delpropensity; j++) {
    for (int i = 1-delpropensity; i <= nx_local+delpropensity; i++) {
      buftmp[m++] = mask[i][j];
    }
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;

      if (screen) {
	m = 0;
	nxtmp = buftmp[m++];
	nytmp = buftmp[m++];
	nztmp = buftmp[m++];
	nxhtmp = buftmp[m++];
	nyhtmp = buftmp[m++];
	nzhtmp = buftmp[m++];
	nxotmp = buftmp[m++];
	nyotmp = buftmp[m++];
	nzotmp = buftmp[m++];
	fprintf(screen,"iproc = %d \n",iproc);
	fprintf(screen,"nxlocal = %d \nnylocal = %d \n",nxtmp,nytmp);
	fprintf(screen,"nx_offset = %d \nny_offset = %d \n",nxotmp,nyotmp);
	m = nrecv;
	for (int j = nytmp+delpropensity; j >= 1-delpropensity; j--) {
	  m-=nxtmp+2*delpropensity;
	  for (int i = 1-delpropensity; i <= nxtmp+delpropensity; i++) {
	    fprintf(screen,"%3d",buftmp[m++]);
	  }
	  fprintf(screen,"\n");
	  m-=nxtmp+2*delpropensity;
	}
      }
      
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
}

/* ----------------------------------------------------------------------
   this is to prevent clustering for undefined child apps
   should eventually replace with pure virtual function
------------------------------------------------------------------------- */

void AppLattice2d::push_connected_neighbors(int i, int j, int** cluster_ids,
					    int id, std::stack<int>*)
{
  error->all("Connectivity not defined for this AppLattice child class");
}


/* ----------------------------------------------------------------------
   this is to prevent clustering for undefined child apps
   should eventually replace with pure virtual function
------------------------------------------------------------------------- */

void AppLattice2d::connected_ghosts(int i, int j, int** cluster_ids,
				    Cluster* clustlist, int idoffset)
{
  error->all("Connectivity not defined for this AppLattice child class");
}
