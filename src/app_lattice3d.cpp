/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppLattice3d::AppLattice3d(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  ntimestep = 0;
  time = 0.0;
  stats_delta = 0.0;
  dump_delta = 0.0;
  temperature = 0.0;
  maxdumpbuf = 0;
  dumpbuf = NULL;
  fp = NULL;
  propensity = NULL;
  site2ijk = NULL;
  ijk2site = NULL;

  // These are sensible defaults, but each app-specific values
  // should be provided in the child constructor.
  dellocal = 0;
  delghost = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice3d::~AppLattice3d()
{
  delete [] dumpbuf;
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

  // app-specific initialization

  init_app();

  // error checks

  if (sweep == NULL && solve == NULL)
    error->all("Lattice app needs a solver or sweeper");

  if (sweep && ((SweepLattice3d *) sweep)->Lkmc && solve == NULL)
    error->all("Must define solver with KMC sweeper");

  if (solve && sweep && ((SweepLattice3d *) sweep)->Lkmc == false)
    error->all("Cannot use solver with non-KMC sweeper");

  if (solve && sweep == NULL && nprocs > 1)
    error->all("Cannot use solver in parallel");

  // initialize arrays
  // propensity only needed if no sweeper
  // if KMC sweep, sweeper will allocate own propensity and site2ijk

  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ijk);
  memory->destroy_3d_T_array(ijk2site);

  int nsites = nx_local*ny_local*nz_local;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "applattice:propensity");
    memory->create_2d_T_array(site2ijk,nsites,3,
			      "applattice2d:site2ijk");
  } else {
    propensity = NULL;
    site2ijk = NULL;
  }

  memory->create_3d_T_array(ijk2site,nx_local+1,ny_local+1,nz_local+1,
			    "applattice3d:ijk2site");

  // initialize lattice <-> site mapping arrays
  // KMC sweeper will overwrite ij2site values

  for (i = 1 ; i <= nx_local; i++)
    for (j = 1 ; j <= ny_local; j++)
      for (k = 1 ; k <= nz_local; k++)
	ijk2site[i][j][k] = (i-1)*ny_local*nz_local + (j-1)*nz_local + k-1;

  if (site2ijk) {
    for (m = 0; m < nsites; m++) {
      i = m / ny_local/nz_local + 1;
      j = (m / nz_local) % ny_local + 1;
      k = m % nz_local + 1;
      site2ijk[m][0] = i;
      site2ijk[m][1] = j;
      site2ijk[m][2] = k;
    }
  }
	 
  // initialize sweeper
  
  if (sweep) sweep->init();

  // initialize propensities for solver
  // if KMC sweep, sweeper does its own init of its propensity arrays

  if (propensity) {
    comm->all(lattice);

    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++)
	for (k = 1 ; k <= nz_local; k++)
	  propensity[ijk2site[i][j][k]] = site_propensity(i,j,k,0);

    solve->init(nsites,propensity);
  }

  // setup future stat and dump calls

  stats_time = time + stats_delta;
  if (stats_delta == 0.0) stats_time = stoptime;
  dump_time = time + dump_delta;
  if (dump_delta == 0.0) dump_time = stoptime;

  // print dump file header and 1st snapshot

  if (dump_delta > 0.0) {
    dump_header();
    dump();
  }

  // print stats header and initial stats
  
  if (me == 0) {
    if (screen) {
      fprintf(screen,"Timestep Time Energy");
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"Timestep Time Energy");
      fprintf(logfile,"\n");
    }
  }

  stats();
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::input_app(char *command, int narg, char **arg)
{
  error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice3d::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

  // init classes used by this app
  
  init();
  timer->init();
  
  // perform the run

  iterate();

  // final statistics
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on solver
 ------------------------------------------------------------------------- */

void AppLattice3d::iterate()
{
  int i,j,k,isite;
  double dt;

  if (time >= stoptime) return

  timer->barrier_start(TIME_LOOP);
  
  int done = 0;
  while (!done) {
    ntimestep++;

    if (propensity) {
      timer->stamp();
      isite = solve->event(&dt);
      timer->stamp(TIME_SOLVE);
      
      if (isite < 0) done = 1;
      else {
	i = site2ijk[isite][0];
	j = site2ijk[isite][1];
	k = site2ijk[isite][2];
	site_event(i,j,k,1);
	time += dt;
	timer->stamp(TIME_APP);
      }

    } else {
      sweep->do_sweep(dt);
      time += dt;
    }

    if (time >= stoptime) done = 1;

    if (time >= stats_time || done) {
      stats();
      stats_time += stats_delta;
      timer->stamp(TIME_OUTPUT);
    }

    if (time >= dump_time || done) {
      if (dump_delta > 0.0) dump();
      dump_time += dump_delta;
      timer->stamp(TIME_OUTPUT);
    }
  }
  
  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice3d::stats()
{
  int i,j,k;
  double energy,all;
  
  comm->all(lattice);

  energy = 0.0;
  for (i = 1; i <= nx_local; i++)
    for (j = 1; j <= ny_local; j++)
      for (k = 1; k <= nz_local; k++)
	energy += site_energy(i,j,k);

  MPI_Allreduce(&energy,&all,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,"%d %f %f\n",ntimestep,time,all);
    if (logfile)
      fprintf(logfile,"%d %f %f\n",ntimestep,time,all);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes element connectivity
   one-time info for viz purposes
------------------------------------------------------------------------- */

void AppLattice3d::dump_header()
{
  // setup comm buf for dumping snapshots

  delete [] dumpbuf;
  maxdumpbuf = 0;
  int mybuf = 2*nx_local*ny_local*nz_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);
  dumpbuf = new int[maxdumpbuf];

  // proc 0 does one-time write of nodes and element connectivity

  if (me) return;

  // number nodes: fast in x, middle in y, slow in z

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestep);
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
  fprintf(fp,"%d\n",ntimestep);
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
------------------------------------------------------------------------- */

void AppLattice3d::dump()
{
  int size_one = 2;

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ELEMENT VALUES\n");
    fprintf(fp,"%d\n",nx_global*ny_global*nz_global);
    fprintf(fp,"ITEM: ELEMENT VALUES\n");
  }

  // pack my lattice values into buffer
  // n = global grid cell (0:Nglobal-1)
  // number cell: fast in x, middle in y, slow in z
  // two triangles per grid cell

  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++)
      for (int k = 1; k <= nz_local; k++) {
	n = (nz_offset+k-1)*ny_global*nx_global + 
	  (ny_offset+j-1)*nx_global + (nx_offset+i-1);
	dumpbuf[m++] = n + 1;
	dumpbuf[m++] = lattice[i][j][k];
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
	MPI_Irecv(dumpbuf,maxdumpbuf,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d\n",dumpbuf[m],dumpbuf[m+1]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dumpbuf,me_size,MPI_INT,0,0,world);
  }
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
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppLattice3d::set_dump(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  if (me == 0) {
    if (fp) fclose(fp);
    fp = fopen(arg[1],"w");
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

  nyz_local = ny_local*nz_local;

  if (dellocal > delghost) 
    error->all("dellocal > delghost: This twisted app could be handled, but not yet");
  if (nx_local < 2*delghost || ny_local < 2*delghost ||
      nz_local < 2*delghost)
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

  nxlo = 1-delghost;
  nxhi = nx_local+delghost;
  nylo = 1-delghost;
  nyhi = ny_local+delghost;
  nzlo = 1-delghost;
  nzhi = nz_local+delghost;
}
