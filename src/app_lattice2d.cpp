/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppLattice2d::AppLattice2d(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
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
}

/* ---------------------------------------------------------------------- */

AppLattice2d::~AppLattice2d()
{
  delete [] dumpbuf;
  memory->sfree(propensity);

  if (fp) {
    fclose(fp);
    fp = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::init()
{
  // error check on other classes

  if (sweep && solve)
    error->all("Lattice app cannot use solver and sweeper");
  if (sweep == NULL && solve == NULL)
    error->all("Lattice app needs a solver or sweeper");

  if (solve && nprocs > 1) error->all("Solvers cannot yet run in parallel");

  // initialize solver:
  //   set propensity of each local site
  //   pass propensity array to solver

  if (solve) {
    if (propensity == NULL) 
      propensity = (double*) memory->smalloc(nx_local*ny_local*sizeof(double),
					     "applattice:propensity");
    comm->all(lattice);

    int i,j,isite;
    for (i = 1 ; i <= nx_local; i++)
      for (j = 1 ; j <= ny_local; j++) {
	isite = ij2site(i,j);
	propensity[isite] = site_propensity(i,j);
      }
    solve->init(nx_local*ny_local,propensity);
  }

  // initialize sweeper

  if (sweep) 
    ((SweepLattice2d *) sweep)->init(nx_global,ny_global,
				     nx_local,ny_local,nx_offset,ny_offset,
				     procwest,proceast,procsouth,procnorth,
				     lattice,temperature);

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

void AppLattice2d::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::input_app(char *command, int narg, char **arg)
{
  error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice2d::run(int narg, char **arg)
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

void AppLattice2d::iterate()
{
  int i,j,k,isite;
  double dt;

  timer->barrier_start(TIME_LOOP);
  
  int done = 0;
  while (!done) {
    ntimestep++;

    if (solve) {
      timer->stamp();
      isite = solve->event(&dt);
      timer->stamp(TIME_SOLVE);
      
      if (isite < 0) done = 1;
      else {
	site2ij(isite,i,j);
	site_event(i,j);
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

void AppLattice2d::stats()
{
  int i,j;
  double energy,all;
  
  comm->all(lattice);

  energy = 0.0;
  for (i = 1; i <= nx_local; i++)
    for (j = 1; j <= ny_local; j++)
      energy += site_energy(i,j);

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

void AppLattice2d::dump_header()
{
  // setup comm buf for dumping snapshots

  delete [] dumpbuf;
  maxdumpbuf = 0;
  int mybuf = 2*nx_local*ny_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);
  dumpbuf = new int[maxdumpbuf];

  // proc 0 does one-time write of nodes and element connectivity

  if (me) return;

  // number nodes: fast in x, slow in y

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestep);
  fprintf(fp,"ITEM: NUMBER OF NODES\n");
  fprintf(fp,"%d\n",(nx_global+1)*(ny_global+1));
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",0.0,(double) nx_global);
  fprintf(fp,"%g %g\n",0.0,(double) ny_global);
  fprintf(fp,"%g %g\n",0.0,0.0);
  fprintf(fp,"ITEM: NODES\n");

  int i,j;
  int m = 0;
  for (j = 0; j <= ny_global; j++)
    for (i = 0; i <= nx_global; i++) {
      m++;
      fprintf(fp,"%d %d %d %d %d\n",m,1,i,j,0);
    }

  // number squares: fast in x, slow in y
  // v1,v2,v3,v4 = 4 corner pts of grid cell in counter-clockwise dir

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestep);
  fprintf(fp,"ITEM: NUMBER OF SQUARES\n");
  fprintf(fp,"%d\n",nx_global*ny_global);
  fprintf(fp,"ITEM: SQUARES\n");

  int v1,v2,v3,v4;
  m = 0;
  for (j = 0; j < ny_global; j++)
    for (i = 0; i < nx_global; i++) {
      v1 = j*(nx_global+1) + i + 1;
      v2 = j*(nx_global+1) + i+1 + 1;
      v3 = (j+1)*(nx_global+1) + i+1 + 1;
      v4 = (j+1)*(nx_global+1) + i + 1;
      m++;
      fprintf(fp,"%d %d %d %d %d %d\n",m,1,v1,v2,v3,v4);
    }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values
------------------------------------------------------------------------- */

void AppLattice2d::dump()
{
  int size_one = 2;

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ELEMENT VALUES\n");
    fprintf(fp,"%d\n",nx_global*ny_global);
    fprintf(fp,"ITEM: ELEMENT VALUES\n");
  }

  // pack my lattice values into buffer
  // n = global grid cell (0:Nglobal-1)
  // number cell: fast in x, slow in y
  // one value per grid cell

  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++) {
      n = (ny_offset+j-1)*nx_global + (nx_offset+i-1);
      dumpbuf[m++] = n + 1;
      dumpbuf[m++] = lattice[i][j];
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
  stats_delta = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::set_dump(int narg, char **arg)
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
   assign nprocs to 2d global lattice so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void AppLattice2d::procs2lattice()
{
  int ipx,ipy,nremain;
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

  if (nx_local < 2 || ny_local < 2)
    error->one("Lattice per proc is too small");

  if (iprocy == 0) procsouth = me + ny_procs - 1;
  else procsouth = me - 1;
  if (iprocy == ny_procs-1) procnorth = me - ny_procs + 1;
  else procnorth = me + 1;

  if (iprocx == 0) procwest = me + nprocs - ny_procs;
  else procwest = me - ny_procs;
  if (iprocx == nx_procs-1) proceast = me - nprocs + ny_procs;
  else proceast = me + ny_procs;
}
