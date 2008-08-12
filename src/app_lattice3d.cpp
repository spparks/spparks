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
#include "app_lattice3d.h"
#include "sweep_lattice3d.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "output.h"

using namespace SPPARKS_NS;

enum {LATFILE,COORDFILE};

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

AppLattice3d::AppLattice3d(SPPARKS *spk, int narg, char **arg) : 
  App(spk,narg,arg)
{
  appclass = LATTICE3D;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  ntimestep = 0;
  time = 0.0;
  temperature = 0.0;
  maxdumpbuf = 0;
  ibufread = NULL;
  ibufdump = NULL;
  dbufdump = NULL;
  fp = NULL;
  propensity = NULL;
  site2ijk = NULL;
  ijk2site = NULL;

  // setup communicator for ghost sites

  comm = new CommLattice3d(spk);

  // app can override these values in its constructor

  delpropensity = 1;
  delevent = 0;
  numrandom = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice3d::~AppLattice3d()
{
  delete comm;

  delete [] ibufread;
  delete [] ibufdump;
  delete [] dbufdump;
  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ijk);
  memory->destroy_3d_T_array(ijk2site);

  if (fp) {
    fclose(fp);
    fp = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::init()
{
  int i,j,k,m;

  // error checks

  if (sweep && strcmp(sweep->style,"lattice3d") != 0)
    error->all("Mismatched sweeper with app lattice");

  if (sweep == NULL && solve == NULL)
    error->all("Lattice app needs a solver or sweeper");

  if (sweep && ((SweepLattice3d *) sweep)->Lkmc && solve == NULL)
    error->all("Must define solver with KMC sweeper");

  if (solve && sweep && ((SweepLattice3d *) sweep)->Lkmc == false)
    error->all("Cannot use solver with non-KMC sweeper");

  if (solve && sweep == NULL && nprocs > 1)
    error->all("Cannot use solver in parallel");

  // app-specific initialization

  init_app();

  // comm init

  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup,
	     delpropensity,delevent);

  // if no sweeper, initialize 3 arrays: propensity, site2i,i2site
  // sweeper allocates its own per-sector versions of these

  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ijk);
  memory->destroy_3d_T_array(ijk2site);

  int nsites = nx_local*ny_local*nz_local;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "app3d:propensity");
    memory->create_2d_T_array(site2ijk,nsites,3,
			      "app3d:site2ijk");
    memory->create_3d_T_array(ijk2site,nx_local+1,ny_local+1,nz_local+1,
			      "app3d:ijk2site");

    for (m = 0; m < nsites; m++) {
      i = m / ny_local/nz_local + 1;
      j = (m / nz_local) % ny_local + 1;
      k = m % nz_local + 1;
      site2ijk[m][0] = i;
      site2ijk[m][1] = j;
      site2ijk[m][2] = k;
    }

    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++)
	for (k = 1 ; k <= nz_local; k++)
	  ijk2site[i][j][k] = (i-1)*ny_local*nz_local + (j-1)*nz_local + k-1;

  } else {
    propensity = NULL;
    site2ijk = NULL;
    ijk2site = NULL;
  }

  // initialize sweeper
  
  if (sweep) {
    sweep->init();
    Lmask = sweep->Lmask;
    mask = ((SweepLattice3d *) sweep)->mask;
  }

  // initialize propensities for KMC solver
  // if KMC sweep, sweeper does its own init of its propensity arrays
  // comm insures ghost sites are set

  if (sweep == NULL) {
    comm->all(lattice);
    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++)
	for (k = 1 ; k <= nz_local; k++)
	  propensity[ijk2site[i][j][k]] = site_propensity(i,j,k);
    solve->init(nsites,propensity);
  }

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) output->set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::input_app(char *command, int narg, char **arg)
{
  error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice3d::run(int narg, char **arg)
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

void AppLattice3d::iterate()
{
  int i,j,k,isite;
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
	i = site2ijk[isite][0];
	j = site2ijk[isite][1];
	k = site2ijk[isite][2];
	site_event(i,j,k,1,random);
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
   update all ghost images of site i,j,k
   called when performing serial KMC on entire domain
------------------------------------------------------------------------- */

void AppLattice3d::update_ghost_sites(int i, int j, int k)
{
  int ijk = lattice[i][j][k];

  // i = 1 plane becomes i = nx_local+1 plane w/ j,k ghosts
  // i = nx_local plane becomes i = 0 plane w/ j,k ghosts

  if (i == 1) {
    lattice[nx_local+1][j][k] = ijk;
    if (j == 1) {
      lattice[nx_local+1][ny_local+1][k] = ijk;
      if (k == 1) lattice[nx_local+1][ny_local+1][nz_local+1] = ijk;
      if (k == nz_local) lattice[nx_local+1][ny_local+1][0] = ijk;
    }
    if (j == ny_local) {
      lattice[nx_local+1][0][k] = ijk;
      if (k == 1) lattice[nx_local+1][0][nz_local+1] = ijk;
      if (k == nz_local) lattice[nx_local+1][0][0] = ijk;
    }
    if (k == 1) lattice[nx_local+1][j][nz_local+1] = ijk;
    if (k == nz_local) lattice[nx_local+1][j][0] = ijk;
  }

  if (i == nx_local) {
    lattice[0][j][k] = ijk;
    if (j == 1) {
      lattice[0][ny_local+1][k] = ijk;
      if (k == 1) lattice[0][ny_local+1][nz_local+1] = ijk;
      if (k == nz_local) lattice[0][ny_local+1][0] = ijk;
    }
    if (j == ny_local) {
      lattice[0][0][k] = ijk;
      if (k == 1) lattice[0][0][nz_local+1] = ijk;
      if (k == nz_local) lattice[0][0][0] = ijk;
    }
    if (k == 1) lattice[0][j][nz_local+1] = ijk;
    if (k == nz_local) lattice[0][j][0] = ijk;
  }

  // j = 1 plane becomes j = ny_local+1 plane w/out i ghosts, w/ k ghosts
  // j = ny_local plane becomes j = 0 plane w/out i ghosts, w/ k ghosts

  if (j == 1) {
    lattice[i][ny_local+1][k] = ijk;
    if (k == 1) lattice[i][ny_local+1][nz_local+1] = ijk;
    if (k == nz_local) lattice[i][ny_local+1][0] = ijk;
  }
  if (j == ny_local) {
    lattice[i][0][k] = ijk;
    if (k == 1) lattice[i][0][nz_local+1] = ijk;
    if (k == nz_local) lattice[i][0][0] = ijk;
  }

  // k = 1 plane becomes k = nz_local+1 plane w/out i,j ghosts
  // k = nz_local plane becomes k = 0 plane w/out i,j ghosts

  if (k == 1) lattice[i][j][nz_local+1] = ijk;
  if (k == ny_local) lattice[i][j][0] = ijk;
}

/* ----------------------------------------------------------------------
   add i,j,k value to site list if different than oldstate and existing sites
------------------------------------------------------------------------- */

void AppLattice3d::add_unique(int oldstate, int &nevent, int *sites,
			      int i, int j, int k)
{
  int value = lattice[i][j][k];
  if (value == oldstate) return;
  for (int m = 0; m < nevent; m++)
    if (value == sites[m]) return;
  sites[nevent++] = value;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice3d::stats(char *strtmp)
{
  int ntimestepall;
  MPI_Allreduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,world);
  sprintf(strtmp," %10d %10g",ntimestepall,time);
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice3d::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s","Step","Time");
}

/* ----------------------------------------------------------------------
   proc 0 writes element connectivity
   one-time info for viz purposes
------------------------------------------------------------------------- */

void AppLattice3d::dump_header()
{
  // setup comm buf for dumping snapshots

  maxdumpbuf = 0;
  int mybuf = nx_local*ny_local*nz_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);

  if (dump_style == LATFILE) {
    delete [] ibufdump;
    ibufdump = new int[2*maxdumpbuf];
  } else {
    delete [] dbufdump;
    dbufdump = new double[5*maxdumpbuf];
  }

  // no header info in file if style = COORDFILE

  if (dump_style == COORDFILE) return;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 does one-time write of nodes and element connectivity

  if (me) return;

  // number nodes: fast in x, middle in y, slow in z

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestepall);
  fprintf(fp,"ITEM: NUMBER OF NODES\n");
  fprintf(fp,"%d\n",(nx_global+1)*(ny_global+1)*(nz_global+1));
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",0.0,(double) nx_global);
  fprintf(fp,"%g %g\n",0.0,(double) ny_global);
  fprintf(fp,"%g %g\n",0.0,(double) nz_global);
  fprintf(fp,"ITEM: NODES\n");

  int i,j,k;
  int m = 0;
  for (k = 0; k <= nz_global; k++)
    for (j = 0; j <= ny_global; j++)
      for (i = 0; i <= nx_global; i++) {
	m++;
	fprintf(fp,"%d %d %d %d %d\n",m,1,i,j,k);
      }

  // v1,v2,v3,v4,v5,v6,v7,v8 = 8 corner pts of grid cell
  // v1-4 are lower plane in counter-clockwise dir, v5-8 are upper plane

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestepall);
  fprintf(fp,"ITEM: NUMBER OF CUBES\n");
  fprintf(fp,"%d\n",nx_global*ny_global*nz_global);
  fprintf(fp,"ITEM: CUBES\n");

  int v1,v2,v3,v4,v5,v6,v7,v8;
  m = 0;
  for (k = 0; k < nz_global; k++)
    for (j = 0; j < ny_global; j++)
      for (i = 0; i < nx_global; i++) {
	v1 = k*(ny_global+1)*(nx_global+1) + j*(nx_global+1) + i + 1;
	v2 = k*(ny_global+1)*(nx_global+1) + j*(nx_global+1) + i+1 + 1;
	v3 = k*(ny_global+1)*(nx_global+1) + (j+1)*(nx_global+1) + i+1 + 1;
	v4 = k*(ny_global+1)*(nx_global+1) + (j+1)*(nx_global+1) + i + 1;
	v5 = (k+1)*(ny_global+1)*(nx_global+1) + j*(nx_global+1) + i + 1;
	v6 = (k+1)*(ny_global+1)*(nx_global+1) + j*(nx_global+1) + i+1 + 1;
	v7 = (k+1)*(ny_global+1)*(nx_global+1) + (j+1)*(nx_global+1) + i+1 + 1;
	v8 = (k+1)*(ny_global+1)*(nx_global+1) + (j+1)*(nx_global+1) + i + 1;
	m++;
	fprintf(fp,"%d %d %d %d %d %d %d %d %d %d\n",
		m,1,v1,v2,v3,v4,v5,v6,v7,v8);
      }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values
   either as lattice sites or as atom coords
------------------------------------------------------------------------- */

void AppLattice3d::dump()
{
  if (dump_style == LATFILE) dump_lattice();
  else dump_coord();
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values, one ELEMENT per site
------------------------------------------------------------------------- */

void AppLattice3d::dump_lattice()
{
  int size_one = 2;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestepall);
    fprintf(fp,"ITEM: NUMBER OF ELEMENT VALUES\n");
    fprintf(fp,"%d\n",nx_global*ny_global*nz_global);
    fprintf(fp,"ITEM: ELEMENT VALUES\n");
  }

  // pack my lattice values into buffer
  // n = global grid cell (0:Nglobal-1)

  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++)
      for (int k = 1; k <= nz_local; k++) {
	n = (nz_offset+k-1)*ny_global*nx_global + 
	  (ny_offset+j-1)*nx_global + (nx_offset+i-1);
	ibufdump[m++] = n + 1;
	ibufdump[m++] = lattice[i][j][k];
      }
  int me_size = m;

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(ibufdump,size_one*maxdumpbuf,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d\n",ibufdump[m],ibufdump[m+1]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(ibufdump,me_size,MPI_INT,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice coords values, one ATOM per site
   app can provide the coords
------------------------------------------------------------------------- */

void AppLattice3d::dump_coord()
{
  int size_one = 5;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestepall);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nx_global*ny_global*nz_global);

    double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;
    box_bounds(&boxxlo,&boxxhi,&boxylo,&boxyhi,&boxzlo,&boxzhi);

    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxxlo,boxxhi);
    fprintf(fp,"%g %g\n",boxylo,boxyhi);
    fprintf(fp,"%g %g\n",boxzlo,boxzhi);
    fprintf(fp,"ITEM: ATOMS\n");
  }

  // pack my lattice coords into buffer
  // n = global grid cell (0:Nglobal-1)

  double x,y,z;
  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++)
      for (int k = 1; k <= nz_local; k++) {
	n = (nz_offset+k-1)*ny_global*nx_global + 
	  (ny_offset+j-1)*nx_global + (nx_offset+i-1);
	dbufdump[m++] = n + 1;
	dbufdump[m++] = lattice[i][j][k];
	xyz(i,j,k,&x,&y,&z);
	dbufdump[m++] = x;
	dbufdump[m++] = y;
	dbufdump[m++] = z;
      }
  int me_size = m;

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(dbufdump,size_one*maxdumpbuf,MPI_DOUBLE,iproc,0,world,
		  &request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d %g %g %g\n",
		static_cast<int> (dbufdump[m]),
		static_cast<int> (dbufdump[m+1]),
		dbufdump[m+2],dbufdump[m+3],dbufdump[m+4]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbufdump,me_size,MPI_DOUBLE,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   box bounds defaults to lattice size
------------------------------------------------------------------------- */

void AppLattice3d::box_bounds(double *xlo, double *xhi, double *ylo,
			      double *yhi, double *zlo, double *zhi)
{
  *xlo = 0.0;
  *ylo = 0.0;
  *zlo = 0.0;
  *xhi = nx_global;
  *yhi = ny_global;
  *zhi = nz_global;
}

/* ----------------------------------------------------------------------
   convert local i,j,k to x,y,z
------------------------------------------------------------------------- */

void AppLattice3d::xyz(int i, int j, int k, double *x, double *y, double *z)
{
  *x = i-1 + nx_offset;
  *y = j-1 + ny_offset;
  *z = k-1 + nz_offset;
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::set_stats(int narg, char **arg)
{
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"your_option_here") == 0) {
      iarg++;
      if (iarg < narg) {
	int itmp = atoi(arg[iarg]);
      } else {
	error->all("Illegal stats command");
      }
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::set_dump(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal dump command");

  if (strcmp(arg[1],"lattice") == 0) dump_style = LATFILE;
  else if (strcmp(arg[1],"coord") == 0) dump_style = COORDFILE;
  else error->all("Illegal dump command");

  if (me == 0) {
    if (fp) fclose(fp);
    fp = fopen(arg[2],"w");
    if (!fp) error->one("Cannot open dump file");
  }
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d global lattice so as to minimize surf area per proc
------------------------------------------------------------------------- */

void AppLattice3d::procs2lattice()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;
  double bestsurf = 2 * (nx_global*ny_global + ny_global*nz_global + 
			 nz_global*nx_global);
  
  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      nremain = nprocs/ipx;
      ipy = 1;
      while (ipy <= nremain) {
        if (nremain % ipy == 0) {
          ipz = nremain/ipy;
	  boxx = float(nx_global)/ipx;
	  boxy = float(ny_global)/ipy;
	  boxz = float(nz_global)/ipz;
	  surf = boxx*boxy + boxy*boxz + boxz*boxx;
	  if (surf < bestsurf) {
	    bestsurf = surf;
	    nx_procs = ipx;
	    ny_procs = ipy;
	    nz_procs = ipz;
	  }
	}
	ipy++;
      }
    }
    ipx++;
  }

  int nyz_procs = ny_procs*nz_procs;
  int iprocx = (me/nyz_procs) % nx_procs;
  nx_offset = iprocx*nx_global/nx_procs;
  nx_local = (iprocx+1)*nx_global/nx_procs - nx_offset;

  int iprocy = (me/nz_procs) % ny_procs;
  ny_offset = iprocy*ny_global/ny_procs;
  ny_local = (iprocy+1)*ny_global/ny_procs - ny_offset;

  int iprocz = (me/1) % nz_procs;
  nz_offset = iprocz*nz_global/nz_procs;
  nz_local = (iprocz+1)*nz_global/nz_procs - nz_offset;

  if (delevent > delpropensity) 
    error->all("Delevent > delpropensity");
  if (nx_local < 2*delpropensity || ny_local < 2*delpropensity ||
      nz_local < 2*delpropensity)
    error->one("Lattice per proc is too small");

  if (iprocz == 0) procdown = me + nz_procs - 1;
  else procdown = me - 1;
  if (iprocz == nz_procs-1) procup = me - nz_procs + 1;
  else procup = me + 1;
    
  if (iprocy == 0) procsouth = me + nyz_procs - nz_procs;
  else procsouth = me - nz_procs;
  if (iprocy == ny_procs-1) procnorth = me - nyz_procs + nz_procs;
  else procnorth = me + nz_procs;
    
  if (iprocx == 0) procwest = me + nprocs - nyz_procs;
  else procwest = me - nyz_procs;
  if (iprocx == nx_procs-1) proceast = me - nprocs + nyz_procs;
  else proceast = me + nyz_procs;

  nxlo = 1-delpropensity;
  nxhi = nx_local+delpropensity;
  nylo = 1-delpropensity;
  nyhi = ny_local+delpropensity;
  nzlo = 1-delpropensity;
  nzhi = nz_local+delpropensity;
}

/* ----------------------------------------------------------------------
   read lattice spin values from file
   unique grain identifier and spin values (agg3d pic file)
   the grain identifier is ignored, only the spin value is stored.
 ------------------------------------------------------------------------- */

void AppLattice3d::read_spins(const char* infile)
{
  char line[MAXLINE];
  char *eof;
  int i,j,k,ii,jj,kk,isite,nbuf,nglobal,ndata,maxbuf,ierr;
  int* ipnt;
  int i1,i2;
  FILE* fpspins;

  // open file

  if (me == 0) {
    fpspins = fopen(infile,"r");
    if (fpspins == NULL) {
      char str[128];
      sprintf(str,"Cannot open file %s",infile);
      error->one(str);
    }
  }

  nglobal = nx_global*ny_global*nz_global;
  ndata = 1;
  maxbuf = 1000;
  delete [] ibufread;
  ibufread = new int[ndata*maxbuf];
  isite = 0;
  while (isite < nglobal) {
    nbuf = 0;
    ipnt = ibufread;
    while (nbuf < maxbuf && isite < nglobal) {
      if (me == 0) {
	eof = fgets(line,MAXLINE,fpspins);
	if (eof == NULL) error->one("Unexpected end of lattice spin file");
	while (line[0]=='#') {
	  eof = fgets(line,MAXLINE,fpspins);
	  if (eof == NULL) error->one("Unexpected end of lattice spin file");
	}

	// Ignore first field if two fields present
	ierr = sscanf(line,"%d %d",&i1,&i2);
	if (ierr == 1) {
	  *ipnt = i1;
	} else if (ierr == 2) {
	  *ipnt = i2;
	} else {
	  error->one("Unexpected value in spin file");
	}
      }
      ipnt+=ndata;
      nbuf++;
      isite++;
    }
    MPI_Bcast(ibufread,ndata*nbuf,MPI_INT,0,world);

    ipnt = ibufread;
    for (int m = isite-nbuf; m < isite; m++) {
      // This file format breaks the SPParKS rule 
      // In agg3d, the innermost loop is on the first coordinate
      i = m % nx_global + 1;
      j = (m / nx_global) % ny_global + 1;
      k = m / ny_global/nx_global + 1;
      ii = i - nx_offset;
      jj = j - ny_offset;
      kk = k - nz_offset;
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local &&
	  kk >= 1 && kk <= nz_local)
	lattice[ii][jj][kk] = ipnt[0];
      ipnt+=ndata;
    }
  }

  if (me == 0) {
    fclose(fpspins);
  }
}

/* ----------------------------------------------------------------------
   this is to prevent clustering for undefined child apps
   should eventually replace with pure virtual function
------------------------------------------------------------------------- */

void AppLattice3d::push_connected_neighbors(int i, int j, int k,
					    int*** cluster_ids,
					    int id, std::stack<int>*)
{
  error->all("Connectivity not defined for this AppLattice child class");
}

/* ----------------------------------------------------------------------
   this is to prevent clustering for undefined child apps
   should eventually replace with pure virtual function
------------------------------------------------------------------------- */

void AppLattice3d::connected_ghosts(int i, int j, int k,
				    int*** cluster_ids,
				    Cluster* clustlist, int idoffset)
{
  error->all("Connectivity not defined for this AppLattice child class");
}
