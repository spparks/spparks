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

enum DumpStyles {LATTICE,COORD};

#define MAXLINE 256

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

  dellocal = 0;
  delghost = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice3d::~AppLattice3d()
{

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

  // app-specific initialization

  init_app();

  // comm init

  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup,
	     delghost,dellocal);

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

  // initialize arrays
  // propensity only needed if no sweeper
  // KMC sweeper will allocate own propensity and site2ijk arrays

  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ijk);
  memory->destroy_3d_T_array(ijk2site);

  int nsites = nx_local*ny_local*nz_local;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "applattice3d:propensity");
    memory->create_2d_T_array(site2ijk,nsites,3,
			      "applattice3d:site2ijk");
  } else {
    propensity = NULL;
    site2ijk = NULL;
  }

  memory->create_3d_T_array(ijk2site,nx_local+1,ny_local+1,nz_local+1,
			    "applattice3d:ijk2site");

  // initialize lattice <-> site mapping arrays
  // they map proc's entire 3d sub-domain to 1d sites and vice versa
  // KMC sweeper will overwrite app's ijk2site values
  // KMC sweeper will create sector-specific site2ijk values and ignore app's

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

  stats_header();
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
  double ptot;
  
  comm->all(lattice);

  energy = 0.0;
  for (i = 1; i <= nx_local; i++)
    for (j = 1; j <= ny_local; j++)
      for (k = 1; k <= nz_local; k++)
	energy += site_energy(i,j,k);

  MPI_Allreduce(&energy,&all,1,MPI_DOUBLE,MPI_SUM,world);

  if (solve == NULL) {
    ptot = 0.0;
  } else {
    ptot = solve->get_total_propensity();
  }

  if (me == 0) {
    if (screen)
      fprintf(screen,"%d %f %f %f\n",ntimestep,time,all,ptot);
    if (logfile)
      fprintf(logfile,"%d %f %f %f\n",ntimestep,time,all,ptot);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice3d::stats_header()
{

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Timestep Time Energy Propensity");
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"Timestep Time Energy Propensity");
      fprintf(logfile,"\n");
    }
  }
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

  if (dump_style == LATTICE) {
    delete [] ibufdump;
    ibufdump = new int[2*maxdumpbuf];
  } else {
    delete [] dbufdump;
    dbufdump = new double[5*maxdumpbuf];
  }

  // no header info in file if style = COORD

  if (dump_style == COORD) return;

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
   either as lattice sites or as atom coords
------------------------------------------------------------------------- */

void AppLattice3d::dump()
{
  if (dump_style == LATTICE) dump_lattice();
  else dump_coord();
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values, one ELEMENT per site
------------------------------------------------------------------------- */

void AppLattice3d::dump_lattice()
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

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
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
	MPI_Irecv(dbufdump,size_one*maxdumpbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d %g %g %g\n",
		static_cast<int> (dbufdump[m]),static_cast<int> (dbufdump[m+1]),
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
  *y = j;
  *z = k;
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
  if (narg != 3) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  if (strcmp(arg[1],"lattice") == 0) dump_style = LATTICE;
  else if (strcmp(arg[1],"coord") == 0) dump_style = COORD;
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

  nx_sector_lo = 1;
  nx_sector_hi = nx_local;
  ny_sector_lo = 1;
  ny_sector_hi = ny_local;
  nz_sector_lo = 1;
  nz_sector_hi = nz_local;
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
  // open file

  if (me == 0) {
    fp = fopen(infile,"r");
    if (fp == NULL) {
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
	eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of lattice spin file");
	while (line[0]=='#') {
	  eof = fgets(line,MAXLINE,fp);
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
    fclose(fp);
  }

}
