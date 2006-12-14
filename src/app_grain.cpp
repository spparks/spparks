/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_grain.h"
#include "comm_grain.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;
using namespace std;

/* ---------------------------------------------------------------------- */

AppGrain::AppGrain(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
{
  int i,j,ii,jj;
  
  if (narg != 5) error->all("Invalid app_style grain command");
  
  dimension = 2;
  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  seed = atoi(arg[4]);
  
  // define proc layout on lattice with 4 periodic neighbors
  
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  
  procs2lattice();
  
  if (me % ny_procs == 0) procsouth = me + ny_procs - 1;
  else procsouth = me - 1;
  if (me % ny_procs == ny_procs-1) procnorth = me - ny_procs + 1;
  else procnorth = me + 1;
  if (me/ny_procs == 0) procwest = me + nprocs - ny_procs;
  else procwest = me - ny_procs;
  if (me/ny_procs == nx_procs-1) proceast = me - nprocs + ny_procs;
  else proceast = me + ny_procs;
  
  // partition lattice across procs
  // allocate my sub-section
  
  int iprocx = me/ny_procs;
  nx_offset = iprocx*nx_global/nx_procs;
  nx_local = (iprocx+1)*nx_global/nx_procs - nx_offset;
  int iprocy = me % ny_procs;
  ny_offset = iprocy*ny_global/ny_procs;
  ny_local = (iprocy+1)*ny_global/ny_procs - ny_offset;
  
  if (nx_local < 2 || ny_local < 2)
    error->one("Lattice per proc is too small");
  
  nx_half = nx_local/2 + 1;
  ny_half = ny_local/2 + 1;
  
  lattice = memory->create_2d_int_array(nx_local+2,ny_local+2,
					"grain:lattice");

  // init ghost spins (take this out later)

  for (i = 0; i <= nx_local+1; i++) 
    for (j = 0; j <= ny_local+1; j++) 
      lattice[i][j] = 0;
  
  // initialize spins
  
  random = new RandomPark(seed);
  for (i = 0; i < 100; i++) random->uniform();
  
  for (i = 1; i <= nx_local; i++) 
    for (j = 1; j <= ny_local; j++) 
      lattice[i][j] = random->irandom(nspins);

  // setup other classes
  
  comm = new CommGrain(spk);
  comm->setup(nx_local,ny_local,nx_half,ny_half,
	      procwest,proceast,procsouth,procnorth);
  
  // default settings

  ntimestep = 0;
  nstats = ndump = 0;
  temperature = 0.0;
  maxbuf = 0;
  buf = NULL;
}

/* ---------------------------------------------------------------------- */

AppGrain::~AppGrain()
{
  memory->destroy_2d_int_array(lattice);
  delete random;
  delete comm;
}

/* ---------------------------------------------------------------------- */

void AppGrain::init()
{
  // sweeping extents for each quadrant
  
  quad[0].xlo = 1;
  quad[0].xhi = nx_half-1;
  quad[0].ylo = 1;
  quad[0].yhi = ny_half-1;

  quad[1].xlo = 1;
  quad[1].xhi = nx_half-1;
  quad[1].ylo = ny_half;
  quad[1].yhi = ny_local;

  quad[2].xlo = nx_half;
  quad[2].xhi = nx_local;
  quad[2].ylo = 1;
  quad[2].yhi = ny_half-1;
  
  quad[3].xlo = nx_half;
  quad[3].xhi = nx_local;
  quad[3].ylo = ny_half;
  quad[3].yhi = ny_local;
 

  // Print layout info
  
  if (me == 0) {
    if (screen) {
      fprintf(screen," nx_global = %d \n ny_global = %d \n",
            nx_global,ny_global);
      fprintf(screen," nx_procs = %d \n ny_procs = %d \n",
            nx_procs,ny_procs);
    }
  }
 
  // setup future stat and dump calls

  stats_next = nstats;
  dump_next = ndump;

  // print dump file header and 1st snapshot

  if (ndump) dump_header();
  if (ndump) dump();

  // print stats header
  
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

  // Print initial stats and dump
  if (nstats) {
    stats();
  }
}

/* ---------------------------------------------------------------------- */

void AppGrain::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);
  else error->all("Invalid command");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppGrain::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  nsweep = atoi(arg[0]);
  
  // error check
  
  //if (solve == NULL) error->all("No solver class defined");
  
  // init classes used by this app
  
  int i;
  double *tmp;
  
  init();
  timer->init();
  
  // perform the run
  
  iterate();
  
  // final statistics
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on sweep solver
 ------------------------------------------------------------------------- */

void AppGrain::iterate()
{
  timer->barrier_start(TIME_LOOP);
  
  for (int i = 0; i < nsweep; i++) {
    ntimestep++;

    /*
    printf("AAA %d %d\n",nx_global,ny_global);
    for (int m = 0; m < nx_global+2; m++)
      printf("BBB %d %d %d %d %d %d %d %d %d %d %d %d\n",
	     lattice[m][0],lattice[m][1],lattice[m][2],lattice[m][3],
	     lattice[m][4],lattice[m][5],lattice[m][6],lattice[m][7],
	     lattice[m][8],lattice[m][9],lattice[m][10],lattice[m][11]);
    */

    timer->stamp();
    comm->communicate(lattice,0);
    timer->stamp(TIME_COMM);

    /*
    printf("AAA %d %d\n",nx_global,ny_global);
    for (int m = 0; m < nx_global+2; m++)
      printf("BBB %d %d %d %d %d %d %d %d %d %d %d %d\n",
	     lattice[m][0],lattice[m][1],lattice[m][2],lattice[m][3],
	     lattice[m][4],lattice[m][5],lattice[m][6],lattice[m][7],
	     lattice[m][8],lattice[m][9],lattice[m][10],lattice[m][11]);
    */

    timer->stamp();
    sweep(0);
    timer->stamp(TIME_SOLVE);
    
    timer->stamp();
       comm->communicate(lattice,1);
    timer->stamp(TIME_COMM);
    
    timer->stamp();
    sweep(1);
    timer->stamp(TIME_SOLVE);
    
    timer->stamp();
       comm->communicate(lattice,2);
    timer->stamp(TIME_COMM);
    
    timer->stamp();
    sweep(2);
    timer->stamp(TIME_SOLVE);
    
    timer->stamp();
    comm->communicate(lattice,3);
    timer->stamp(TIME_COMM);
    
    timer->stamp();
    sweep(3);
    timer->stamp(TIME_SOLVE);
    
    if (ntimestep == stats_next) {
      timer->stamp();
      stats();
      stats_next += nstats;
      timer->stamp(TIME_OUTPUT);
    }
    if (ntimestep == dump_next) {
      timer->stamp();
      dump();
      dump_next += ndump;
      timer->stamp(TIME_OUTPUT);
    }
  }
  
  timer->barrier_stop(TIME_LOOP);
}

/* ---------------------------------------------------------------------- */
   
void AppGrain::sweep(int index)
{
  int i,j,iold,inew,nold,nnew;

  int xlo = quad[index].xlo;
  int xhi = quad[index].xhi;
  int ylo = quad[index].ylo;
  int yhi = quad[index].yhi;

  for (i = xlo; i <= xhi; i++) {
    for (j = ylo; j <= yhi; j++) {
      iold = lattice[i][j];
      nold = 0;
      if (iold == lattice[i-1][j-1]) nold++;
      if (iold == lattice[i-1][j]) nold++;
      if (iold == lattice[i-1][j+1]) nold++;
      if (iold == lattice[i][j-1]) nold++;
      if (iold == lattice[i][j+1]) nold++;
      if (iold == lattice[i+1][j-1]) nold++;
      if (iold == lattice[i+1][j]) nold++;
      if (iold == lattice[i+1][j+1]) nold++;

      inew = random->irandom(nspins);
      nnew = 0;
      if (inew == lattice[i-1][j-1]) nnew++;
      if (inew == lattice[i-1][j]) nnew++;
      if (inew == lattice[i-1][j+1]) nnew++;
      if (inew == lattice[i][j-1]) nnew++;
      if (inew == lattice[i][j+1]) nnew++;
      if (inew == lattice[i+1][j-1]) nnew++;
      if (inew == lattice[i+1][j]) nnew++;
      if (inew == lattice[i+1][j+1]) nnew++;

      if (nold <= nnew) lattice[i][j] = inew;
    }
  }
}

/* ----------------------------------------------------------------------
   Compute total energy of system, quadrant by quandrant
------------------------------------------------------------------------- */

double AppGrain::compute_energy()
{
  double energy_local,energy_global;
  energy_local = 0.0;

  for (int index = 0; index < 4; index++) {

//     // Dump spins before communicate
//     char *fstring = "Quadrant %4d, before communicate, in compute_energy()";
//     int len = strlen(fstring)+32;
//     char *title = new char[len];
//     sprintf(title,fstring,index);
//     dump_detailed(title);
//     delete [] title;

    comm->communicate(lattice,index);

//     // Dump spins after communicate
//     fstring = "Quadrant %4d, after communicate, in compute_energy()";
//     len = strlen(fstring)+32;
//     title = new char[len];
//     sprintf(title,fstring,index);
//     dump_detailed(title);
//     delete [] title;

    energy_local += energy_quadrant(index);
  }

  MPI_Allreduce(&energy_local,&energy_global,1,MPI_DOUBLE,MPI_SUM,world);
    
  return energy_global;
}

/* ----------------------------------------------------------------------
   Update all site in one quadrant
------------------------------------------------------------------------- */

double AppGrain::energy_quadrant(int index)
{
  int i,j,ik,nk;
  double energy;

  energy = 0.0;

  for (i = quad[index].xlo; i <= quad[index].xhi; i++) {
    for (j = quad[index].ylo; j <= quad[index].yhi; j++) {
      ik = lattice[i][j];
      nk = 0;
      if (ik != lattice[i-1][j-1]) nk++;
      if (ik != lattice[i-1][j]) nk++;
      if (ik != lattice[i-1][j+1]) nk++;
      if (ik != lattice[i][j-1]) nk++;
      if (ik != lattice[i][j+1]) nk++;
      if (ik != lattice[i+1][j-1]) nk++;
      if (ik != lattice[i+1][j]) nk++;
      if (ik != lattice[i+1][j+1]) nk++;
      energy+=nk;
    }
  }

  return energy;

}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppGrain::stats()
{
  double energy;

  energy = compute_energy();
  if (me == 0) {
    if (screen) {
      fprintf(screen,"%d %f \n",ntimestep,energy);
    }
    if (logfile) {
      fprintf(logfile,"%d %f \n",ntimestep,energy);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes element connectivity
   one-time info for viz purposes
------------------------------------------------------------------------- */

void AppGrain::dump_header()
{
  // setup comm buf for dumping snapshots

  delete [] buf;
  maxbuf = 0;
  int mybuf = 4*nx_local*ny_local;
  MPI_Allreduce(&mybuf,&maxbuf,1,MPI_INT,MPI_MAX,world);
  buf = new int[maxbuf];
  
  // proc 0 does one-time write of nodes and element connectivity

  if (me) return;

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
  for (i = 0; i <= nx_global; i++)
    for (j = 0; j <= ny_global; j++) {
      m++;
      fprintf(fp,"%d %d %d %d %d\n",m,1,i,j,0);
    }

  // v1,v2,v3,v4 = 4 corner pts of grid cell, in counter-clockwise dir

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestep);
  fprintf(fp,"ITEM: NUMBER OF TRIANGLES\n");
  fprintf(fp,"%d\n",2*nx_global*ny_global);
  fprintf(fp,"ITEM: TRIANGLES\n");

  int v1,v2,v3,v4;
  m = 0;
  for (i = 0; i < nx_global; i++)
    for (j = 0; j < ny_global; j++) {
      v1 = i*(ny_global+1) + j + 1;
      v2 = (i+1)*(ny_global+1) + j + 1;
      v3 = (i+1)*(ny_global+1) + j+1 + 1;
      v4 = i*(ny_global+1) + j+1 + 1;
      m++;
      fprintf(fp,"%d %d %d %d %d\n",m,1,v1,v3,v2);
      m++;
      fprintf(fp,"%d %d %d %d %d\n",m,1,v1,v4,v3);
    }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice spin values
------------------------------------------------------------------------- */

void AppGrain::dump()
{
  int size_one = 2;

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ELEMENT VALUES\n");
    fprintf(fp,"%d\n",2*nx_global*ny_global);
    fprintf(fp,"ITEM: ELEMENT VALUES\n");
  }

  // pack my lattice values into buffer
  // n = global grid cell (0:Nglobal-1)
  // two triangles per grid cell

  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++) {
      n = (nx_offset+i-1)*ny_global + (ny_offset+j-1);
      buf[m++] = 2*n + 1;
      buf[m++] = lattice[i][j];
      buf[m++] = 2*n + 2;
      buf[m++] = lattice[i][j];
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
	MPI_Irecv(buf,maxbuf,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fp,"%d %d\n",buf[m],buf[m+1]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,me_size,MPI_INT,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of grains to the screen in 2D layout
   all the spins for each processor domain are printed out
------------------------------------------------------------------------- */

void AppGrain::dump_detailed(char* title)
{
  int nsend,nrecv,nxtmp,nytmp,nxhtmp,nyhtmp,nxotmp,nyotmp;
  int size_one = 1;

  // set up communication buffer
  // maxbuf must equal the maximum number of spins on one domain 
  // plus some extra stuff
  maxbuf = ((nx_global-1)/nx_procs+3)*((ny_global-1)/ny_procs+3)+6;
  buf = (int*) memory->smalloc(maxbuf*sizeof(int),"appgrain:dump_detailed:buf");

  nsend = (nx_local+2)*(ny_local+2)+6;
  if (maxbuf < nsend) {
    error->one("maxbuf size too small in AppGrain::dump_detailed()");
  }

  // proc 0 writes interactive dump header

  if (me == 0) {
    fprintf(screen,"*** Interactive Dump ***\n");
    fprintf(screen,"Title = %s\n",title);
    fprintf(screen,"nx_global = %d ny_global = %d \n",nx_global,ny_global);
  }

  int m = 0;

  // pack local layout info into buffer

  buf[m++] = nx_local;
  buf[m++] = ny_local;
  buf[m++] = nx_half;
  buf[m++] = ny_half;
  buf[m++] = nx_offset;
  buf[m++] = ny_offset;

  // pack my lattice values into buffer
  // Need to violate normal ordering in order to simplify output
  for (int j = 0; j <= ny_local+1; j++) {
    for (int i = 0; i <= nx_local+1; i++) {
      buf[m++] = lattice[i][j];
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
	MPI_Irecv(buf,maxbuf,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;

      m = 0;
      nxtmp = buf[m++];
      nytmp = buf[m++];
      nxhtmp = buf[m++];
      nyhtmp = buf[m++];
      nxotmp = buf[m++];
      nyotmp = buf[m++];
      fprintf(screen,"iproc = %d \n",iproc);
      fprintf(screen,"nxlocal = %d \nnylocal = %d \n",nxtmp,nytmp);
      fprintf(screen,"nx_half = %d \nny_half = %d \n",nxhtmp,nyhtmp);
      fprintf(screen,"nx_offset = %d \nny_offset = %d \n",nxotmp,nyotmp);
      m = nrecv;
      for (int j = nytmp+1; j >= 0; j--) {
	m-=nxtmp+2;
	for (int i = 0; i <= nxtmp+1; i++) {
	  fprintf(screen,"%3d",buf[m++]);
	}
	fprintf(screen,"\n");
	m-=nxtmp+2;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buf);
}

/* ---------------------------------------------------------------------- */

void AppGrain::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppGrain::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  nstats = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppGrain::set_dump(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal dump command");
  ndump = atoi(arg[0]);
  if (me == 0) {
    fp = fopen(arg[1],"w");
    if (fp == NULL) error->one("Cannot open dump file");
  }
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d nx/nyglobal lattice so as to minimize perimeter 
------------------------------------------------------------------------- */

void AppGrain::procs2lattice()
{
  int ipx,ipy,nremain;
  // Need float arithmetic to handle uneven ratios
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
  
}
