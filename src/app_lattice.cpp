/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_lattice.h"
#include "sweep_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

enum{NONE,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       RANDOM_2D,RANDOM_3D,FILENAME};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppLattice::AppLattice(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
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
  dbuf = NULL;
  fp = NULL;
  propensity = NULL;

  // app can override these values in its constructor

  dellocal = 0;
  delghost = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  delete [] dbuf;
  memory->sfree(propensity);
  delete [] latfile;

  memory->sfree(id);
  memory->destroy_2d_T_array(xyz);
  memory->sfree(numneigh);
  memory->destroy_2d_T_array(neighbor);

  if (fp) {
    fclose(fp);
    fp = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::options(int narg, char **arg)
{
  // defaults

  latstyle = NONE;
  latfile = NULL;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lattice") == 0) {
      if (iarg+1 > narg) error->all("Illegal app_style command");
      if (strcmp(arg[iarg+1],"sq/4n") == 0) latstyle = SQ_4N;
      else if (strcmp(arg[iarg+1],"sq/8n") == 0) latstyle = SQ_8N;
      else if (strcmp(arg[iarg+1],"tri") == 0) latstyle = TRI;
      else if (strcmp(arg[iarg+1],"sc/6n") == 0) latstyle = SC_6N;
      else if (strcmp(arg[iarg+1],"sc/26n") == 0) latstyle = SC_26N;
      else if (strcmp(arg[iarg+1],"fcc") == 0) latstyle = FCC;
      else if (strcmp(arg[iarg+1],"bcc") == 0) latstyle = BCC;
      else if (strcmp(arg[iarg+1],"diamond") == 0) latstyle = DIAMOND;
      else if (strcmp(arg[iarg+1],"random/2d") == 0) latstyle = RANDOM_2D;
      else if (strcmp(arg[iarg+1],"random/3d") == 0) latstyle = RANDOM_3D;
      else if (strcmp(arg[iarg+1],"file") == 0) latstyle = FILENAME;
      else error->all("Illegal app_style command");

      if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
	if (iarg+5 > narg) error->all("Illegal app_style command");
	dimension = 2;
	latconst = atof(arg[iarg+2]);
	nx = atoi(arg[iarg+3]);
	ny = atoi(arg[iarg+4]);
	iarg += 5;
      } else if (latstyle == SC_6N || latstyle == SC_26N ||
		 latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
	if (iarg+6 > narg) error->all("Illegal app_style command");
	dimension = 3;
	latconst = atof(arg[iarg+2]);
	nx = atoi(arg[iarg+3]);
	ny = atoi(arg[iarg+4]);
	nz = atoi(arg[iarg+5]);
	iarg += 6;
      } else if (latstyle == RANDOM_2D) {
	if (iarg+6 > narg) error->all("Illegal app_style command");
	dimension = 2;
	nrandom = atoi(arg[iarg+2]);
	xprd = atof(arg[iarg+3]);
	yprd = atof(arg[iarg+4]);
	cutoff = atof(arg[iarg+5]);
	iarg += 6;
      } else if (latstyle == RANDOM_3D) {
	if (iarg+7 > narg) error->all("Illegal app_style command");
	dimension = 3;
	nrandom = atoi(arg[iarg+2]);
	xprd = atof(arg[iarg+3]);
	yprd = atof(arg[iarg+4]);
	zprd = atof(arg[iarg+4]);
	cutoff = atof(arg[iarg+6]);
	iarg += 7;
      } else if (latstyle == FILENAME) {
	if (iarg+3 > narg) error->all("Illegal app_style command");
	int n = strlen(arg[2]) + 1;
	latfile = new char[n];
	latfile = strcpy(latfile,arg[2]);
	iarg += 3;
      }
    } else error->all("Illegal app_style command ");
  }

  if (latstyle == NONE) error->all("Illegal app_style command");
}

/* ----------------------------------------------------------------------
   generate lattice and processor decomposition
 ------------------------------------------------------------------------- */

void AppLattice::create_lattice()
{
  if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI || 
      latstyle == SC_6N || latstyle == SC_26N || 
      latstyle == FCC || latstyle == BCC || latstyle == DIAMOND)
    structured_lattice();
  else if (latstyle == RANDOM_2D || latstyle == RANDOM_3D)
    random_lattice();
  else if (latstyle == FILENAME)
    file_lattice();
}

/* ----------------------------------------------------------------------
   structured lattice
 ------------------------------------------------------------------------- */

void AppLattice::structured_lattice()
{
  // determine box extent

  if (latstyle == SQ_4N || latstyle == SQ_8N) {
    boxxlo = boxylo = 0.0;
    boxxhi = nx * latconst;
    boxyhi = ny * latconst;
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == TRI) {
    boxxlo = boxylo = 0.0;
    boxxhi = nx * latconst;
    boxyhi = ny * latconst * sqrt(3.0);
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == SC_6N || latstyle == SC_26N || 
	     latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
    boxxlo = boxylo = boxzlo = 0.0;
    boxxhi = nx * latconst;
    boxyhi = ny * latconst;
    boxzhi = nz * latconst;
  }

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  // partition domain

  if (dimension == 2) procs2lattice_2d();
  else if (dimension == 3) procs2lattice_3d();

  // basis atoms of each unit cell depend on lattice

  double basis[8][3];

  if (latstyle == SQ_4N || latstyle == SQ_8N) {
    nbasis = 1;
    nglobal = nbasis * nx*ny;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
  } else if (latstyle == TRI) {
    nbasis = 2;
    nglobal = nbasis * nx*ny;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
    basis[1][0] = 0.5; basis[1][1] = 0.5; basis[1][2] = 0.0;

  } else if (latstyle == SC_6N || latstyle == SC_26N) {
    nbasis = 1;
    nglobal = nbasis * nx*ny*nz;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
  } else if (latstyle == FCC) {
    nbasis = 4;
    nglobal = nbasis * nx*ny*nz;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
    basis[1][0] = 0.0; basis[1][1] = 0.5; basis[1][2] = 0.5;
    basis[2][0] = 0.5; basis[2][1] = 0.0; basis[2][2] = 0.5;
    basis[3][0] = 0.5; basis[3][1] = 0.5; basis[3][2] = 0.0;
  } else if (latstyle == BCC) {
    nbasis = 2;
    nglobal = nbasis * nx*ny*nz;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
    basis[1][0] = 0.5; basis[1][1] = 0.5; basis[1][2] = 0.5;
  } else if (latstyle == DIAMOND) {
    nbasis = 8;
    nglobal = nbasis * nx*ny*nz;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
    basis[1][0] = 0.0; basis[1][1] = 0.5; basis[1][2] = 0.5;
    basis[2][0] = 0.5; basis[2][1] = 0.0; basis[2][2] = 0.5;
    basis[3][0] = 0.5; basis[3][1] = 0.5; basis[3][2] = 0.0;
    basis[4][0] = 0.25; basis[4][1] = 0.25; basis[4][2] = 0.25;
    basis[5][0] = 0.25; basis[5][1] = 0.75; basis[5][2] = 0.75;
    basis[6][0] = 0.75; basis[6][1] = 0.25; basis[6][2] = 0.75;
    basis[7][0] = 0.75; basis[7][1] = 0.75; basis[7][2] = 0.25;
  }

  // generate lattice of sites
  // 1st pass = count lattice points I own in my sub-domain
  // 2nd pass = generate xyz coords and store them with site ID

  if (dimension == 2) nz = 1;

  double x,y,z;
  double latconstx = latconst;
  double latconsty = latconst;
  double latconstz = latconst;
  if (latstyle == TRI) latconsty *= sqrt(3.0);

  nlocal = 0;
  int i,j,k,m;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
	for (m = 0; m < nbasis; m++) {
	  x = (i + basis[m][0]) * latconstx;
	  y = (j + basis[m][1]) * latconsty;
	  if (dimension == 2) z = 0.0;
	  else z = (k + basis[m][2]) * latconstz;
	  if (x < subxlo || x >= subxhi || 
	      y < subylo || y >= subyhi || 
	      z < subzlo || z >= subzhi) continue;
	  nlocal++;
	}

  id = (int *) memory->smalloc(nlocal*sizeof(int),"app:id");
  memory->create_2d_T_array(xyz,nlocal,3,"app:xyz");

  nlocal = 0;
  int count = 0;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
	for (m = 0; m < nbasis; m++) {
	  count++;
	  x = (i + basis[m][0]) * latconstx;
	  y = (j + basis[m][1]) * latconsty;
	  if (dimension == 2) z = 0.0;
	  else z = (k + basis[m][2]) * latconstz;

	  if (x < subxlo || x >= subxhi || 
	      y < subylo || y >= subyhi || 
	      z < subzlo || z >= subzhi) continue;

	  id[nlocal] = count;
	  xyz[nlocal][0] = x;
	  xyz[nlocal][1] = y;
	  xyz[nlocal][2] = z;
	  nlocal++;
	}

  // allocate arrays to store lattice connectivity

  if (latstyle == SQ_4N) maxconnect = 4;
  else if (latstyle == SQ_8N) maxconnect = 8;
  else if (latstyle == TRI) maxconnect = 6;
  else if (latstyle == SC_6N) maxconnect = 6;
  else if (latstyle == SC_26N) maxconnect = 26;
  else if (latstyle == FCC) maxconnect = 12;
  else if (latstyle == BCC) maxconnect = 8;
  else if (latstyle == DIAMOND) maxconnect = 4;

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  memory->create_2d_T_array(neighbor,nlocal,maxconnect,"app:neighbor");

  // create connectivity offsets

  memory->create_3d_T_array(cmap,nbasis,maxconnect,4,"app:cmap");
  offsets();

  // generate lattice connectivity for each site
  // connect() computes global index of Jth neighbor of global site I
  // global index = 1 to Nglobal

  for (i = 0; i < nlocal; i++) {
    numneigh[i] = maxconnect;
    for (j = 0; j < maxconnect; j++) {
      neighbor[i][j] = connect(id[i],j);
      if (neighbor[i][j] <= 0 || neighbor[i][j] > nglobal)
	error->all("Bad connectivity result");

      // connectivity check on distance

      /*
      double dx = xyz[i][0] - xyz[neighbor[i][j]-1][0];
      double dy = xyz[i][1] - xyz[neighbor[i][j]-1][1];
      double dz = xyz[i][2] - xyz[neighbor[i][j]-1][2];
      if (dx > 0.5*xprd) dx -= xprd;
      if (dx < -0.5*xprd) dx += xprd;
      if (dy > 0.5*yprd) dy -= yprd;
      if (dy < -0.5*yprd) dy += yprd;
      if (dz > 0.5*zprd) dz -= zprd;
      if (dz < -0.5*zprd) dz += zprd;
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      printf("DIST %d %d: %g\n",id[i],id[neighbor[i][j]-1],r);
      */
    }
  }

  memory->destroy_3d_T_array(cmap);

  // map global to local indices via maps
  // needs to be a parallel operation that generates ghost sites

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < maxconnect; j++)
      neighbor[i][j]--;

  nghost = 0;
}

/* ----------------------------------------------------------------------
   random lattice
 ------------------------------------------------------------------------- */

void AppLattice::random_lattice()
{
  int i,j,n;

  if (latstyle == RANDOM_2D) {
    boxxlo = boxylo = 0.0;
    boxxhi = xprd;
    boxyhi = yprd;
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == RANDOM_3D) {
    boxxlo = boxylo = boxzlo = 0.0;
    boxxhi = xprd;
    boxyhi = yprd;
    boxzhi = zprd;
  }

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  // partition domain

  if (dimension == 2) procs2lattice_2d();
  else if (dimension == 3) procs2lattice_3d();

  nglobal = nrandom;

  // generate random sites
  // 1st pass = count sites I own in my sub-domain
  // 2nd pass = generate xyz coords and store them with site ID
  // save and restore seed between passes

  double x,y,z;
  int seed = random->seed;

  nlocal = 0;
  for (n = 0; n < nglobal; n++) {
    x = xprd * random->uniform();
    y = yprd * random->uniform();
    if (dimension == 3) z = zprd * random->uniform();
    else z = 0.0;
    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    nlocal++;
  }

  random->seed = seed;

  id = (int *) memory->smalloc(nlocal*sizeof(int),"app:id");
  memory->create_2d_T_array(xyz,nlocal,3,"app:xyz");

  nlocal = 0;
  for (n = 0; n < nglobal; n++) {
    x = xprd * random->uniform();
    y = yprd * random->uniform();
    if (dimension == 3) z = zprd * random->uniform();
    else z = 0.0;

    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    
    id[nlocal] = n + 1;
    xyz[nlocal][0] = x;
    xyz[nlocal][1] = y;
    xyz[nlocal][2] = z;
    nlocal++;
  }

  // count max neighbors thru expensive N^2 loop (SERIAL ONLY)

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  double dx,dy,dz,rsq;
  double cutsq = cutoff*cutoff;

  for (i = 0; i < nlocal; i++)
    for (j = i+1; j < nlocal; j++) {
      dx = xyz[i][0] - xyz[j][0];
      dy = xyz[i][1] - xyz[j][1];
      dz = xyz[i][2] - xyz[j][2];

      if (dx < -0.5*xprd) dx += xprd;
      else if (dx > 0.5*xprd) dx -= xprd;
      if (dy < -0.5*yprd) dy += yprd;
      else if (dy > 0.5*yprd) dy -= yprd;
      if (dimension == 3) {
	if (dz < -0.5*zprd) dz += zprd;
	else if (dz > 0.5*zprd) dz -= zprd;
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutsq) {
	numneigh[i]++;
	numneigh[j]++;
      }
    }

  // allocate arrays to store lattice connectivity

  maxconnect = 0;
  for (i = 0; i < nlocal; i++) maxconnect = MAX(maxconnect,numneigh[i]);
  if (maxconnect == 0) error->all("Random lattice has no connectivity");

  memory->create_2d_T_array(neighbor,nlocal,maxconnect,"app:neighbor");

  // generate lattice connectivity for each site thru expensive N^2 loop
  // SERIAL ONLY

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  for (i = 0; i < nlocal; i++)
    for (j = i+1; j < nlocal; j++) {
      dx = xyz[i][0] - xyz[j][0];
      dy = xyz[i][1] - xyz[j][1];
      dz = xyz[i][2] - xyz[j][2];

      if (dx < -0.5*xprd) dx += xprd;
      else if (dx > 0.5*xprd) dx -= xprd;
      if (dy < -0.5*yprd) dy += yprd;
      else if (dy > 0.5*yprd) dy -= yprd;
      if (dimension == 3) {
	if (dz < -0.5*zprd) dz += zprd;
	else if (dz > 0.5*zprd) dz -= zprd;
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutsq) {
	neighbor[i][numneigh[i]++] = j;
	neighbor[j][numneigh[j]++] = i;
      }
    }

  nghost = 0;
}

/* ----------------------------------------------------------------------
   read lattice from file
 ------------------------------------------------------------------------- */

void AppLattice::file_lattice()
{
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  int i,j,m;

  // app-specific initialization

  init_app();

  // error checks

  if (sweep && strcmp(sweep->style,"lattice") != 0)
    error->all("Mismatched sweeper with app lattice");

  if (sweep == NULL && solve == NULL)
    error->all("Lattice app needs a solver or sweeper");

  if (sweep && ((SweepLattice*) sweep)->Lkmc && solve == NULL)
    error->all("Must define solver with KMC sweeper");

  if (solve && sweep && ((SweepLattice*) sweep)->Lkmc == false)
    error->all("Cannot use solver with non-KMC sweeper");

  if (solve && sweep == NULL && nprocs > 1)
    error->all("Cannot use solver in parallel");

  // initialize arrays
  // propensity only needed if no sweeper
  // if KMC sweep, sweeper will allocate own propensity and site2ij

  memory->sfree(propensity);

  if (sweep == NULL)
    propensity = (double*) memory->smalloc(nlocal*sizeof(double),
					   "app:propensity");
  else propensity = NULL;

  // initialize sweeper
  
  if (sweep) sweep->init();

  // initialize propensities for solver
  // if KMC sweep, sweeper does its own init of its propensity arrays

  if (propensity) {
    comm->all(lattice);

    for (i = 0 ; i < nlocal; i++)
      propensity[i] = site_propensity(i,0);

    solve->init(nlocal,propensity);
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

void AppLattice::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice::input_app(char *command, int narg, char **arg)
{
  error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppLattice::run(int narg, char **arg)
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

void AppLattice::iterate()
{
  int i,j,isite;
  double dt;

  if (time >= stoptime) return;

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
	site_event(isite,1);
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

void AppLattice::stats()
{
  int i,j;
  double energy,all;
  
  comm->all(lattice);

  energy = 0.0;
  for (i = 0; i < nlocal; i++)
    energy += site_energy(i);

  MPI_Allreduce(&energy,&all,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,"%d %f %f\n",ntimestep,time,all);
    if (logfile)
      fprintf(logfile,"%d %f %f\n",ntimestep,time,all);
  }

}

/* ----------------------------------------------------------------------
   allocate dump buffer
------------------------------------------------------------------------- */

void AppLattice::dump_header()
{
  // setup comm buf for dumping snapshots

  delete [] dbuf;
  maxdumpbuf = 0;
  MPI_Allreduce(&nlocal,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);

  dbuf = new double[5*maxdumpbuf];
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values as atom coords
------------------------------------------------------------------------- */

void AppLattice::dump()
{
  int size_one = 5;

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nglobal);
    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxxlo,boxxhi);
    fprintf(fp,"%g %g\n",boxylo,boxyhi);
    fprintf(fp,"%g %g\n",boxzlo,boxzhi);
    fprintf(fp,"ITEM: ATOMS\n");
  }

  // pack my lattice coords into buffer

  double x,y,z;
  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    dbuf[m++] = id[i];
    dbuf[m++] = lattice[i];
    dbuf[m++] = xyz[i][0];
    dbuf[m++] = xyz[i][1];
    dbuf[m++] = xyz[i][2];
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
	MPI_Irecv(dbuf,size_one*maxdumpbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d %g %g %g\n",
		static_cast<int> (dbuf[m]),static_cast<int> (dbuf[m+1]),
		dbuf[m+2],dbuf[m+3],dbuf[m+4]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbuf,me_size,MPI_DOUBLE,0,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_dump(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  if (strcmp(arg[1],"coord") != 0) error->all("Illegal dump command");

  if (me == 0) {
    if (fp) fclose(fp);
    fp = fopen(arg[2],"w");
    if (!fp) error->one("Cannot open dump file");
  }
}

/* ----------------------------------------------------------------------
   assign nprocs to global box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void AppLattice::procs2lattice_2d()
{
  int ipx,ipy,nremain;
  double boxx,boxy,surf;
  double bestsurf = 2.0 * (xprd * yprd);
  
  // loop thru all possible factorizations of nprocs
  // surf = perimeter of a proc sub-domain
 
  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = xprd/ipx;
      boxy = yprd/ipy;
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
  int iprocy = me % ny_procs;

  if (iprocy == 0) procsouth = me + ny_procs - 1;
  else procsouth = me - 1;
  if (iprocy == ny_procs-1) procnorth = me - ny_procs + 1;
  else procnorth = me + 1;

  if (iprocx == 0) procwest = me + nprocs - ny_procs;
  else procwest = me - ny_procs;
  if (iprocx == nx_procs-1) proceast = me - nprocs + ny_procs;
  else proceast = me + ny_procs;

  subxlo = iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = (iprocx+1) * xprd/nx_procs;
  else subxhi = xprd;

  subylo = iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = (iprocy+1) * yprd/ny_procs;
  else subyhi = yprd;

  subzlo = -0.5;
  subzhi = 0.5;
}

/* ---------------------------------------------------------------------- */

void AppLattice::procs2lattice_3d()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;
  double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);
  
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
	  boxx = xprd/ipx;
	  boxy = yprd/ipy;
	  boxz = zprd/ipz;
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

  int nyz_procs = ny_procs * nz_procs;

  int iprocx = (me/(nyz_procs)) % nx_procs;
  int iprocy = (me/nz_procs) % ny_procs;
  int iprocz = (me/1) % nz_procs;

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

  subxlo = iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = (iprocx+1) * xprd/nx_procs;
  else subxhi = xprd;

  subylo = iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = (iprocy+1) * yprd/ny_procs;
  else subyhi = yprd;

  subzlo = iprocz * zprd/nz_procs;
  if (iprocz < nz_procs-1) subzhi = (iprocz+1) * zprd/nz_procs;
  else subzhi = zprd;
}

/* ---------------------------------------------------------------------- */

int AppLattice::connect(int iglobal, int ineigh)
{
  int i,j,k,m,n;

  if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
    i = (iglobal-1)/nbasis % nx;
    j = (iglobal-1)/nbasis / nx;
    m = (iglobal-1) % nbasis;
    i += cmap[m][ineigh][0];
    j += cmap[m][ineigh][1];
    m = cmap[m][ineigh][2];
    if (i < 0) i += nx;
    if (i >= nx) i -= nx;
    if (j < 0) j += ny;
    if (j >= ny) j -= ny;
    n = j*nx*nbasis + i*nbasis + m + 1;

  } else if (latstyle == SC_6N || latstyle == SC_26N || 
	     latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
    i = (iglobal-1)/nbasis % nx;
    j = (iglobal-1)/nbasis / nx % ny;;
    k = (iglobal-1)/nbasis / (nx*ny);
    m = (iglobal-1) % nbasis;
    i += cmap[m][ineigh][0];
    j += cmap[m][ineigh][1];
    k += cmap[m][ineigh][2];
    m = cmap[m][ineigh][3];
    if (i < 0) i += nx;
    if (i >= nx) i -= nx;
    if (j < 0) j += ny;
    if (j >= ny) j -= ny;
    if (k < 0) k += nz;
    if (k >= nz) k -= nz;
    n = k*nx*ny*nbasis + j*nx*nbasis + i*nbasis + m + 1;
  }

  return n;
}

/* ---------------------------------------------------------------------- */

void AppLattice::offsets()
{
  if (latstyle == SQ_4N) {
    cmap[0][0][0] = -1; cmap[0][0][1] =  0; cmap[0][0][2] = 0;
    cmap[0][1][0] =  1; cmap[0][1][1] =  0; cmap[0][1][2] = 0;
    cmap[0][2][0] =  0; cmap[0][2][1] = -1; cmap[0][2][2] = 0;
    cmap[0][3][0] =  0; cmap[0][3][1] =  1; cmap[0][3][2] = 0;

  } else if (latstyle == SQ_8N) {
    cmap[0][0][0] = -1; cmap[0][0][1] =  0; cmap[0][0][2] = 0;
    cmap[0][1][0] =  1; cmap[0][1][1] =  0; cmap[0][1][2] = 0;
    cmap[0][2][0] =  0; cmap[0][2][1] = -1; cmap[0][2][2] = 0;
    cmap[0][3][0] =  0; cmap[0][3][1] =  1; cmap[0][3][2] = 0;
    cmap[0][4][0] = -1; cmap[0][4][1] = -1; cmap[0][4][2] = 0;
    cmap[0][5][0] = -1; cmap[0][5][1] =  1; cmap[0][5][2] = 0;
    cmap[0][6][0] =  1; cmap[0][6][1] = -1; cmap[0][6][2] = 0;
    cmap[0][7][0] =  1; cmap[0][7][1] =  1; cmap[0][7][2] = 0;

  } else if (latstyle == TRI) {
    cmap[0][0][0] = -1; cmap[0][0][1] =  0; cmap[0][0][2] = 0;
    cmap[0][1][0] =  1; cmap[0][1][1] =  0; cmap[0][1][2] = 0;
    cmap[0][2][0] =  0; cmap[0][2][1] =  0; cmap[0][2][2] = 1;
    cmap[0][3][0] = -1; cmap[0][3][1] =  0; cmap[0][3][2] = 1;
    cmap[0][4][0] = -1; cmap[0][4][1] = -1; cmap[0][4][2] = 1;
    cmap[0][5][0] =  0; cmap[0][5][1] = -1; cmap[0][5][2] = 1;

    cmap[0][0][1] = -1; cmap[0][0][1] =  0; cmap[0][0][2] = 1;
    cmap[1][1][0] =  1; cmap[1][1][1] =  0; cmap[1][1][2] = 1;
    cmap[1][2][0] =  0; cmap[1][2][1] =  0; cmap[1][2][2] = 0;
    cmap[1][3][0] =  0; cmap[1][3][1] =  1; cmap[1][3][2] = 0;
    cmap[1][4][0] =  0; cmap[1][4][1] =  1; cmap[1][4][2] = 0;
    cmap[1][5][0] =  1; cmap[1][5][1] =  1; cmap[1][5][2] = 0;

  } else if (latstyle == SC_6N) {
    cmap[0][0][0] = -1; cmap[0][0][1] =  0; cmap[0][0][2] =  0;
    cmap[0][1][0] =  1; cmap[0][1][1] =  0; cmap[0][1][2] =  0;
    cmap[0][2][0] =  0; cmap[0][2][1] = -1; cmap[0][2][2] =  0;
    cmap[0][3][0] =  0; cmap[0][3][1] =  1; cmap[0][3][2] =  0;
    cmap[0][4][0] =  0; cmap[0][4][1] =  0; cmap[0][4][2] =  -1;
    cmap[0][5][0] =  0; cmap[0][5][1] =  0; cmap[0][5][2] =   1;
    for (int m = 0; m < nbasis; m++) cmap[0][m][3] = 0;
 
  } else if (latstyle == SC_26N) {
    int m = 0;
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	for (int k = -1; k <= 1; k++) {
	  if (i == 0 && j == 0 && k == 0) continue;
	  cmap[0][m][0] = i; cmap[0][m][1] = j; cmap[0][m][2] = k;
	  cmap[0][m][3] = 0;
	  m++;
	}

  } else if (latstyle == FCC) {
    int m = 0;                             // 1st basis atom
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same x
      for (int jj = -1; jj <= 0; jj++) {
	cmap[0][m][0] = 0; cmap[0][m][1] = ii; cmap[0][m][2] = jj;
	cmap[0][m][3] = 1;
	m++;
      }
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same y
      for (int jj = -1; jj <= 0; jj++) {
	cmap[0][m][0] = ii; cmap[0][m][1] = 0; cmap[0][m][2] = jj;
	cmap[0][m][3] = 2;
	m++;
      }
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same y
      for (int jj = -1; jj <= 0; jj++) {
	cmap[0][m][0] = ii; cmap[0][m][1] = jj; cmap[0][m][2] = 0;
	cmap[0][m][3] = 3;
	m++;
      }

    m = 0;                                 // 2nd basis atom
    for (int ii = 0; ii <= 1; ii++)        // 4 neighs with same x
      for (int jj = 0; jj <= 1; jj++) {
	cmap[1][m][0] = 0; cmap[1][m][1] = ii; cmap[1][m][2] = jj;
	cmap[1][m][3] = 0;
	m++;
      }
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same y
      for (int jj = 0; jj <= 1; jj++) {
	cmap[1][m][0] = ii; cmap[1][m][1] = 0; cmap[1][m][2] = jj;
	cmap[1][m][3] = 3;
	m++;
      }
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same z
      for (int jj = 0; jj <= 1; jj++) {
	cmap[1][m][0] = ii; cmap[1][m][1] = jj; cmap[1][m][2] = 0;
	cmap[1][m][3] = 2;
	m++;
      }

    m = 0;                                 // 3rd basis atom
    for (int ii = -1; ii <= 0; ii++)       // 4 neighs with same x
      for (int jj = 0; jj <= 1; jj++) {
	cmap[2][m][0] = 0; cmap[2][m][1] = ii; cmap[2][m][2] = jj;
	cmap[2][m][3] = 3;
	m++;
      }
    for (int ii = 0; ii <= 1; ii++)        // 4 neighs with same y
      for (int jj = 0; jj <= 1; jj++) {
	cmap[2][m][0] = ii; cmap[2][m][1] = 0; cmap[2][m][2] = jj;
	cmap[2][m][3] = 0;
	m++;
      }
    for (int ii = 0; ii <= 1; ii++)       // 4 neighs with same z
      for (int jj = -1; jj <= 0; jj++) {
	cmap[2][m][0] = ii; cmap[2][m][1] = jj; cmap[2][m][2] = 0;
	cmap[2][m][3] = 1;
	m++;
      }

    m = 0;                                 // 4th basis atom
    for (int ii = 0; ii <= 1; ii++)        // 4 neighs with same x
      for (int jj = -1; jj <= 0; jj++) {
	cmap[3][m][0] = 0; cmap[3][m][1] = ii; cmap[3][m][2] = jj;
	cmap[3][m][3] = 2;
	m++;
      }
    for (int ii = 0; ii <= 1; ii++)       // 4 neighs with same y
      for (int jj = -1; jj <= 0; jj++) {
	cmap[3][m][0] = ii; cmap[3][m][1] = 0; cmap[3][m][2] = jj;
	cmap[3][m][3] = 1;
	m++;
      }
    for (int ii = 0; ii <= 1; ii++)        // 4 neighs with same z
      for (int jj = 0; jj <= 1; jj++) {
	cmap[3][m][0] = ii; cmap[3][m][1] = jj; cmap[3][m][2] = 0;
	cmap[3][m][3] = 0;
	m++;
      }

  } else if (latstyle == BCC) {
    int m = 0;
    for (int i = -1; i <= 0; i++)
      for (int j = -1; j <= 0; j++)
	for (int k = -1; k <= 0; k++) {
	  cmap[0][m][0] = i; cmap[0][m][1] = j; cmap[0][m][2] = k;
	  cmap[0][m][3] = 1;
	  m++;
	}

    m = 0;
    for (int i = 0; i <= 1; i++)
      for (int j = 0; j <= 1; j++)
	for (int k = 0; k <= 1; k++) {
	  cmap[1][m][0] = i; cmap[1][m][1] = j; cmap[1][m][2] = k;
	  cmap[1][m][3] = 0;
	  m++;
	}

  } else if (latstyle == DIAMOND) {
    int ibasis[8][3] = {{0,0,0}, {0,2,2}, {2,0,2}, {2,2,0}, 
			{1,1,1}, {1,3,3}, {3,1,3}, {3,3,1}};
    int bond[4][3] = {{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1}};

    int n,m,x,y,z,i,j,k,ib;

    for (n = 0; n < nbasis; n++) {
      for (m = 0; m < maxconnect; m++) {
	if (n < 4) {
	  x = ibasis[n][0] + bond[m][0];
	  y = ibasis[n][1] + bond[m][1];
	  z = ibasis[n][2] + bond[m][2];
	} else {
	  x = ibasis[n][0] - bond[m][0];
	  y = ibasis[n][1] - bond[m][1];
	  z = ibasis[n][2] - bond[m][2];
	}

	if (x < 0) {
	  x += 4;
	  i = -1;
	} else if (x > 3) {
	  x -= 4;
	  i = 1;
	} else i = 0;
	if (y < 0) {
	  y += 4;
	  j = -1;
	} else if (y > 3) {
	  y -= 4;
	  j = 1;
	} else j = 0;
	if (z < 0) {
	  z += 4;
	  k = -1;
	} else if (z > 3) {
	  z -= 4;
	  k = 1;
	} else k = 0;

	cmap[n][m][0] = i; cmap[n][m][1] = j; cmap[n][m][2] = k;
	
	for (ib = 0; ib < nbasis; ib++)
	  if (x == ibasis[ib][0] && y == ibasis[ib][1] && 
	      z == ibasis[ib][2]) cmap[n][m][3] = ib;
      }
    }
  }
}
