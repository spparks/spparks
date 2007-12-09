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

enum{LATTICE,COORD};

#define MAXLINE 256

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
  ibuf = NULL;
  dbuf = NULL;
  fp = NULL;
  propensity = NULL;
  site2ij = NULL;
  ij2site = NULL;

  // setup communicator for ghost sites

  comm = new CommLattice2d(spk);

  // app can override these values in its constructor

  dellocal = 0;
  delghost = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice2d::~AppLattice2d()
{
  delete [] ibuf;
  delete [] dbuf;
  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ij);
  memory->destroy_2d_T_array(ij2site);

  if (fp) {
    fclose(fp);
    fp = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice2d::init()
{
  int i,j,m;

  // app-specific initialization

  init_app();

  // comm init

  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth,
	     delghost,dellocal);

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

  // initialize arrays
  // propensity only needed if no sweeper
  // KMC sweeper will allocate own propensity and site2ij arrays

  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ij);
  memory->destroy_2d_T_array(ij2site);

  int nsites = nx_local*ny_local;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "applattice2d:propensity");
    memory->create_2d_T_array(site2ij,nsites,2,
			      "applattice2d:site2ij");
  } else {
    propensity = NULL;
    site2ij = NULL;
  }

  memory->create_2d_T_array(ij2site,nx_local+1,ny_local+1,
			    "applattice2d:ij2site");

  // initialize lattice <-> site mapping arrays
  // they map proc's entire 2d sub-domain to 1d sites and vice versa
  // KMC sweeper will overwrite app's ij2site values
  // KMC sweeper will create sector-specific site2ij values and ignore app's

  for (i = 1 ; i <= nx_local; i++)
    for (j = 1 ; j <= ny_local; j++)
      ij2site[i][j] = (i-1)*ny_local + j-1;

  if (site2ij) {
    for (m = 0; m < nsites; m++) {
      i = m / ny_local + 1;
      j = m % ny_local + 1;
      site2ij[m][0] = i;
      site2ij[m][1] = j;
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
	propensity[ij2site[i][j]] = site_propensity(i,j,0);

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
  
  stats_header();
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
	i = site2ij[isite][0];
	j = site2ij[isite][1];
	site_event(i,j,1);
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
  double ptot;

  comm->all(lattice);

  energy = 0.0;
  for (i = 1; i <= nx_local; i++)
    for (j = 1; j <= ny_local; j++)
      energy += site_energy(i,j);

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

void AppLattice2d::stats_header()
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

void AppLattice2d::dump_header()
{
  // setup comm buf for dumping snapshots

  delete [] ibuf;
  delete [] dbuf;
  maxdumpbuf = 0;
  int mybuf = nx_local*ny_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);

  if (dump_style == LATTICE) ibuf = new int[2*maxdumpbuf];
  else dbuf = new double[5*maxdumpbuf];

  // no header info in file if style = COORD

  if (dump_style == COORD) return;

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
   either as lattice sites or as atom coords
------------------------------------------------------------------------- */

void AppLattice2d::dump()
{
  if (dump_style == LATTICE) dump_lattice();
  else dump_coord();
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values, one ELEMENT per site
------------------------------------------------------------------------- */

void AppLattice2d::dump_lattice()
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
      ibuf[m++] = n + 1;
      ibuf[m++] = lattice[i][j];
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
	MPI_Irecv(ibuf,size_one*maxdumpbuf,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d\n",ibuf[m],ibuf[m+1]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(ibuf,me_size,MPI_INT,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice coords values, one ATOM per site
   app can provide the coords
------------------------------------------------------------------------- */

void AppLattice2d::dump_coord()
{
  int size_one = 5;

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nx_global*ny_global);

    double boxxlo,boxxhi,boxylo,boxyhi;
    box_bounds(&boxxlo,&boxxhi,&boxylo,&boxyhi);

    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxxlo,boxxhi);
    fprintf(fp,"%g %g\n",boxylo,boxyhi);
    fprintf(fp,"%g %g\n",0.0,0.0);
    fprintf(fp,"ITEM: ATOMS\n");
  }

  // pack my lattice coords into buffer
  // n = global grid cell (0:Nglobal-1)

  double x,y,z;
  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++) {
      n = (ny_offset+j-1)*nx_global + (nx_offset+i-1);
      dbuf[m++] = n + 1;
      dbuf[m++] = lattice[i][j];
      xy(i-1+nx_offset,j-1+ny_offset,&x,&y);
      dbuf[m++] = x;
      dbuf[m++] = y;
      dbuf[m++] = 0.0;
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

/* ----------------------------------------------------------------------
   box bounds defaults to lattice size
------------------------------------------------------------------------- */

void AppLattice2d::box_bounds(double *xlo, double *xhi,
			      double *ylo, double *yhi)
{
  *xlo = 0.0;
  *ylo = 0.0;
  *xhi = nx_global;
  *yhi = ny_global;
}

/* ----------------------------------------------------------------------
   lattice coords defaults to i,j,k
------------------------------------------------------------------------- */

void AppLattice2d::xy(int i, int j, double *x, double *y)
{
  *x = i;
  *y = j;
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

  if (dellocal > delghost) 
    error->all("dellocal > delghost: This twisted app could be handled, but not yet");
  if (nx_local < 2*delghost || ny_local < 2*delghost)
    error->one("Lattice per proc is too small");

  if (iprocy == 0) procsouth = me + ny_procs - 1;
  else procsouth = me - 1;
  if (iprocy == ny_procs-1) procnorth = me - ny_procs + 1;
  else procnorth = me + 1;

  if (iprocx == 0) procwest = me + nprocs - ny_procs;
  else procwest = me - ny_procs;
  if (iprocx == nx_procs-1) proceast = me - nprocs + ny_procs;
  else proceast = me + ny_procs;

  nxlo = 1-delghost;
  nxhi = nx_local+delghost;
  nylo = 1-delghost;
  nyhi = ny_local+delghost;

  nx_sector_lo = 1;
  nx_sector_hi = nx_local;
  ny_sector_lo = 1;
  ny_sector_hi = ny_local;
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
  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delghost)*((ny_global-1)/ny_procs+1+2*delghost)+9;
  nsend = (nx_local+2*delghost)*(ny_local+2*delghost)+9;
  if (maxbuftmp < nsend) 
    error->one("maxbuftmp size too small in AppGrain::dump_detailed()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"applattice2d:dump_detailed:buftmp");

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

  for (int j = 1-delghost; j <= ny_local+delghost; j++) {
    for (int i = 1-delghost; i <= nx_local+delghost; i++) {
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
	for (int j = nytmp+delghost; j >= 1-delghost; j--) {
	  m-=nxtmp+2*delghost;
	  for (int i = 1-delghost; i <= nxtmp+delghost; i++) {
	    fprintf(screen,"%3d",buftmp[m++]);
	  }
	  fprintf(screen,"\n");
	  m-=nxtmp+2*delghost;
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
  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delghost)*((ny_global-1)/ny_procs+1+2*delghost)+9;
  nsend = (nx_local+2*delghost)*(ny_local+2*delghost)+9;
  if (maxbuftmp < nsend) 
    error->one("maxbuftmp size too small in AppGrain::dump_detailed_mask()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"applattice2d:dump_detailed_mask:buftmp");

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

  for (int j = 1-delghost; j <= ny_local+delghost; j++) {
    for (int i = 1-delghost; i <= nx_local+delghost; i++) {
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
	for (int j = nytmp+delghost; j >= 1-delghost; j--) {
	  m-=nxtmp+2*delghost;
	  for (int i = 1-delghost; i <= nxtmp+delghost; i++) {
	    fprintf(screen,"%3d",buftmp[m++]);
	  }
	  fprintf(screen,"\n");
	  m-=nxtmp+2*delghost;
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
   read lattice spin values from file
   unique grain identifier and spin values (agg2d pic file)
   the grain identifier is ignored, only the spin value is stored.
 ------------------------------------------------------------------------- */

void AppLattice2d::read_spins(const char* infile)
{
  char line[MAXLINE];
  char *eof;
  int i,j,ii,jj,isite,nbuf,nglobal,ndata,maxbuf,ierr;
  int* ipnt;
  // open file

  if (me == 0) {
    fp = fopen(infile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open file %s",infile);
      error->one(str);
    }
  }

  nglobal = nx_global*ny_global;
  ndata = 2;
  maxbuf = 1000;
  delete [] ibuf;
  ibuf = new int[ndata*maxbuf];
  isite = 0;
  while (isite < nglobal) {
    nbuf = 0;
    ipnt = ibuf;
    while (nbuf < maxbuf && isite < nglobal) {
      if (me == 0) {
	eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of lattice spin file");
	while (line[0]=='#') {
	  eof = fgets(line,MAXLINE,fp);
	  if (eof == NULL) error->one("Unexpected end of lattice spin file");
	}

	ierr = sscanf(line,"%d %d",ipnt,ipnt+1);
	if (ierr != 2) {
	  error->one("Unexpected value in spin file");
	}
      }
      ipnt+=ndata;
      nbuf++;
      isite++;
    }
    MPI_Bcast(ibuf,ndata*nbuf,MPI_INT,0,world);

    ipnt = ibuf;
    for (int m = isite-nbuf; m < isite; m++) {
      // This file format breaks the SPParKS rule 
      // In agg3d, the innermost loop is on the first coordinate
      i = m % nx_global + 1;
      j = m / nx_global + 1;
      ii = i - nx_offset;
      jj = j - ny_offset;
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local)
	lattice[ii][jj] = ipnt[1];
      ipnt+=ndata;
    }
  }

  if (me == 0) {
    fclose(fp);
  }

  delete [] ibuf;

}
