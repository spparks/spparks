/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_grain_nfw.h"
#include "random_park.h"
#include "solve.h"
#include "comm_grain_2d.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;
using namespace std;

/* ---------------------------------------------------------------------- */

AppGrainNfw::AppGrainNfw(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
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

//   // Detailed dump for debugging

//   char *fstring = "Initial state, in AppGrainNfw() \n";
//   int len = strlen(fstring)+32;
//   char *title = new char[len];
//   sprintf(title,fstring);
//   dump_detailed(title);
//   delete [] title;

//   // setup other classes
  
  comm = new CommGrain2D(spk);
  comm->setup(nx_local,ny_local,nx_half,ny_half,
	      procwest,proceast,procsouth,procnorth);

  // default settings

  ntimestep = 0;
  stats_delta = 0.0;
  dump_delta = 0.0;
  temperature = 0.0;
  maxbuf = 0;
  buf = NULL;
  fp = NULL;
  solve = NULL;
  propensity = NULL;
}

/* ---------------------------------------------------------------------- */

AppGrainNfw::~AppGrainNfw()
{
  memory->destroy_2d_int_array(lattice);
  delete random;
  delete comm;
  delete [] propensity;

  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void AppGrainNfw::init()
{

  // initialize propensities for solver

  init_propensity();

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

  stats_time = time + stats_delta;
  if (stats_delta == 0.0) stats_time = stoptime;
  dump_time = time + dump_delta;
  if (dump_delta == 0.0) dump_time = stoptime;

  // print dump file header and 1st snapshot

  dump_header();
  dump();

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
  stats();
}

/* ---------------------------------------------------------------------- */

void AppGrainNfw::input(char *command, int narg, char **arg)
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

void AppGrainNfw::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);
  
  // init classes used by this app
  
  init();

//   // Detailed dump for debugging

//   char *fstring = "Initial state, in run() \n";
//   int len = strlen(fstring)+32;
//   char *title = new char[len];
//   sprintf(title,fstring);
//   dump_detailed(title);
//   delete [] title;

  timer->init();
  
  // perform the run
  
  iterate();
  
//   // Detailed dump for debugging

//   char *fstring2 = "Final state, in run() \n";
//   int len2 = strlen(fstring2)+32;
//   char *title2 = new char[len2];
//   sprintf(title2,fstring2);
//   dump_detailed(title2);
//   delete [] title2;

  // final statistics
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   Compute total energy of system, quadrant by quandrant
------------------------------------------------------------------------- */

double AppGrainNfw::compute_energy()
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

double AppGrainNfw::energy_quadrant(int index)
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

void AppGrainNfw::stats()
{
  double energy;

  energy = compute_energy();
  if (me == 0) {
    if (screen) {
      fprintf(screen,"%d %f %f\n",ntimestep,time,energy);
    }
    if (logfile) {
      fprintf(logfile,"%d %f %f\n",ntimestep,time,energy);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes element connectivity
   one-time info for viz purposes
------------------------------------------------------------------------- */

void AppGrainNfw::dump_header()
{
  // setup comm buf for dumping snapshots

  if (!fp) return;

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

void AppGrainNfw::dump()
{
  if (!fp) return;

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

void AppGrainNfw::dump_detailed(char* title)
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
    error->one("maxbuf size too small in AppGrainNfw::dump_detailed()");
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

void AppGrainNfw::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppGrainNfw::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppGrainNfw::set_dump(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (me == 0) {
    fp = fopen(arg[1],"w");
    if (!fp) error->one("Cannot open dump file");
  }
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d nx/nyglobal lattice so as to minimize perimeter 
------------------------------------------------------------------------- */

void AppGrainNfw::procs2lattice()
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

/* ----------------------------------------------------------------------
   iterate on N-Fold Way solver
------------------------------------------------------------------------- */

void AppGrainNfw::iterate()
{
  double dt;
  int i,j;
  int done = 0;
  int isite;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    ntimestep++;

    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

//     fprintf(screen,"Picked site %d \n",isite);
//     fprintf(screen,"Picked dt %g \n",dt);

    // Check if solver failed to pick an event

    if (isite < 0) {

      done = 1;
      
    } else {

      // Get index for this site

      i = (isite / ny_local) + 1;
      j = isite - (i-1)*ny_local + 1;

      // Pick new spin from neighbor spins
      
      update_spin(i,j);
    
      timer->stamp(TIME_APP);

      // Update propensity for this site and neighbors
      
      update_propensity(i,j);

      // update time by Gillepsie dt

      time += dt;
      if (time >= stoptime) done = 1;

    }

    timer->stamp();
    if (time > stats_time || done) {
      stats();
      stats_time += stats_delta;
      timer->stamp(TIME_OUTPUT);
    }

    if (time > dump_time || done) {
      dump();
      dump_time += dump_delta;
      timer->stamp(TIME_OUTPUT);
    }

//     char *fstring = "Iteration %d, in iterate() \n";
//     int len = strlen(fstring)+32;
//     char *title = new char[len];
//     sprintf(title,fstring,ntimestep);
//     dump_detailed(title);
//     delete [] title;

  }

  timer->barrier_stop(TIME_LOOP);
}

// Pick new spin from neighbor spins

void AppGrainNfw::update_spin(const int& i, const int& j) {
  int ns,spins[9],nspins[9],ik,jk;
  int icand,ncand,candlist[9];

  // Initialize neighbor spin list
  ik = lattice[i][j];
  spins[0] = ik;
  nspins[0] = 0;
  ns = 1;

  // Survey each neighbor

  jk = lattice[i-1][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i-1][j];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i-1][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);

  // Use survey to identify candidates
  ncand = 0;
  for (int is=1;is<ns;is++) {
    if (nspins[is] >= nspins[0]) {
      candlist[ncand] = spins[is];
      ncand++;
    }
  }

  if (ncand <= 0) error->one("Failed to find candidate spin");
  icand = random->irandom(ncand)-1;

  // Update lattice sites and periodic images
  lattice[i][j] = candlist[icand];

  if (i==1) {
    lattice[nx_local+1][j] = lattice[i][j];
    if (j==1)
      lattice[nx_local+1][ny_local+1] = lattice[i][j];
    if (j==ny_local)
      lattice[nx_local+1][0] = lattice[i][j];
  }

  if (i==nx_local) {
    lattice[0][j] = lattice[i][j];
    if (j==1)
      lattice[0][ny_local+1] = lattice[i][j];
    if (j==ny_local)
      lattice[0][0] = lattice[i][j];
  }

  if (j==1)
    lattice[i][ny_local+1] = lattice[i][j];
  if (j==ny_local)
    lattice[i][0] = lattice[i][j];
}
    
// Update propensity for this site and neighbors

void AppGrainNfw::update_propensity(const int& i, const int& j) const {

  int ndepends = 9;
  int depends[9];
  int k,ii,jj,isite;

  timer->stamp();

  k = 0;

  // For now, we will explitly enforce PBCs.

  ii = i > 1 ? i-1 : nx_local;  
  jj = j > 1 ? j-1 : ny_local;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i > 1 ? i-1 : nx_local;  
  jj = j;
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i > 1 ? i-1 : nx_local;  
  jj = j < ny_local ? j+1 : 1;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i;
  jj = j > 1 ? j-1 : ny_local;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i;
  jj = j;
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i;
  jj = j < ny_local ? j+1 : 1;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i < nx_local ? i+1 : 1;  
  jj = j > 1 ? j-1 : ny_local;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i < nx_local ? i+1 : 1;  
  jj = j;
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  ii = i < nx_local ? i+1 : 1;  
  jj = j < ny_local ? j+1 : 1;  
  isite = (ii-1)*ny_local+jj-1;
  depends[k] = isite;
  propensity[isite] = compute_propensity(ii,jj);
  k++;

  timer->stamp(TIME_APP);


  // update propensities of events

//   fprintf(screen,"i = %d j = %d \n",i,j);
//   fprintf(screen,"ii = %d jj = %d \n",ii,jj);
//   fprintf(screen,"ndepends %d \n",ndepends);
//   for (int idepend = 0; idepend<ndepends; idepend++) 
//     fprintf(screen,"depends[i] %d propensity[i] %g \n",depends[idepend],propensity[idepend]);
  

  solve->update(ndepends,depends,propensity);
  timer->stamp(TIME_UPDATE);

}

// Initialize propensities for all sites

void AppGrainNfw::init_propensity() {
  int n,isite;

//   fprintf(screen,"*** Propensity Dump ***\n");
//   fprintf(screen,"Title = %s\n","No title yet");
//   fprintf(screen,"nx_local = %d ny_local = %d \n",nx_local,ny_local);

  // Compute propensities in temporary array and pass to solve
  n = nx_local*ny_local;
  propensity = (double*) memory->smalloc(n*sizeof(double),"appgrainnfw:init_propensity:propensity");

//   fprintf(screen,"*** Propensity Dump ***\n");
//   fprintf(screen,"Title = %s\n","No title yet");
//   fprintf(screen,"nx_local = %d ny_local = %d \n",nx_local,ny_local);
//   fprintf(screen,"n = %d \n",n);

  for (int i = 1 ; i <= nx_local; i++) { 
    for (int j = 1 ; j <= ny_local; j++) {

      isite = (i-1)*ny_local+j-1;
      propensity[isite] = compute_propensity(i,j);
//       fprintf(screen,"i = %d j = %d isite = %d p = %g \n",i,j,isite,propensity[isite]);

    }
  }
    
  // error check
  
  if (solve == NULL) error->all("No solver class defined");
  
  solve->init(n,propensity);

}

// Compute propensities for site (i,j)

double AppGrainNfw::compute_propensity(const int& i, const int& j) const {
  int ns,spins[9],nspins[9],ik,jk;
  int ncand;
  double prop;

  // Initialize neighbor spin list
  ik = lattice[i][j];
  spins[0] = ik;
  nspins[0] = 0;
  ns = 1;

  // Survey each neighbor

  jk = lattice[i-1][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i-1][j];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i-1][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j-1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j];
  survey_neighbor(ik,jk,ns,spins,nspins);
  jk = lattice[i+1][j+1];
  survey_neighbor(ik,jk,ns,spins,nspins);
  
  // Use survey to identify candidates
  
  ncand = 0;
  for (int is=1;is<ns;is++) {
    if (nspins[is] >= nspins[0]) {
      ncand++;
    }
  }

  prop = ncand;

//   fprintf(screen,"ncand %d \n",ncand);
//   fprintf(screen,"prop %g \n",prop);

  return prop;
}

void AppGrainNfw::survey_neighbor(const int& ik, const int& jk, int& ns, int spins[], int nspins[]) const {
  bool Lfound;
  int is;

  // First check if matches current spin
  if (jk == ik) {
    nspins[0]++;
  // Otherwise compare with other existing spins in list
  } else {
    Lfound = false;
    for (is=1;is<ns;is++) {
      if (jk == spins[is]) {
	Lfound = true;
	break;
      }
    }
    // If found, increment counter 
    if (Lfound) {
      nspins[is]++;
    // If not, create new survey entry
    } else {
      spins[ns] = jk;
      nspins[ns] = 1;
      ns++;
    }
  }

}
    

