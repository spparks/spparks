/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_grain.h"
#include "sweep_grain.h"
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

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (narg < 6) error->all("Invalid app_style grain command");
  
  // default/dummy settings
  ntimestep = 0;
  nstats = ndump = 0;
  temperature = 0.0;
  maxdumpbuf = 0;
  dumpbuf = NULL;
  fp = NULL;

  dimension = 0;
  nx_global = 0;
  ny_global = 0;
  nz_global = 0;
  nx_local = 0;
  ny_local = 0;
  nz_local = 0;
  nx_offset = -1;
  ny_offset = -1;
  nz_offset = -1;
  lattice = NULL;
  lat_2d = NULL;
  lat_3d = NULL;
  procwest = -1;
  proceast = -1;
  procsouth = -1;
  procnorth = -1;
  procup = -1;
  procdown = -1;
  nspins = 0;
  temperature = 0.0;

  int iarg=1;
  if (strcmp(arg[iarg],"2d") == 0) {
    iarg++;
    dimension = 2;
    if (narg != 6) error->all("Invalid app_style grain 2d command");
    nx_global = atoi(arg[iarg++]);
    ny_global = atoi(arg[iarg++]);
  } else if (strcmp(arg[iarg],"3d") == 0) {
    iarg++;
    dimension = 3;
    if (narg != 7) error->all("Invalid app_style grain 3d command");
    nx_global = atoi(arg[iarg++]);
    ny_global = atoi(arg[iarg++]);
    nz_global = atoi(arg[iarg++]);
  } else error->all("Illegal dimension specifier in app_style grain command");
  nspins = atoi(arg[iarg++]);
  seed = atoi(arg[iarg++]);

}

/* ---------------------------------------------------------------------- */

AppGrain::~AppGrain()
{
  memory->destroy_3d_T_array(lattice);
  delete random;

  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void AppGrain::init()
{
  int i,j,k,ii,jj,kk;
  
  // define proc layout on lattice with 4 periodic neighbors

  if (dimension == 2) {
  
    // partition lattice across procs
    // allocate my sub-section
  
    procs2lattice_2d();
    int iprocx = me/ny_procs;
    nx_offset = iprocx*nx_global/nx_procs;
    nx_local = (iprocx+1)*nx_global/nx_procs - nx_offset;
    int iprocy = me % ny_procs;
    ny_offset = iprocy*ny_global/ny_procs;
    ny_local = (iprocy+1)*ny_global/ny_procs - ny_offset;
    
    if (nx_local < 2 || ny_local < 2)
      error->one("Lattice per proc is too small");
    
    // Figure out who neighbor processors are    

    if (me % ny_procs == 0) procsouth = me + ny_procs - 1;
    else procsouth = me - 1;
    if (me % ny_procs == ny_procs-1) procnorth = me - ny_procs + 1;
    else procnorth = me + 1;
    if (me/ny_procs == 0) procwest = me + nprocs - ny_procs;
    else procwest = me - ny_procs;
    if (me/ny_procs == nx_procs-1) proceast = me - nprocs + ny_procs;
    else proceast = me + ny_procs;

    memory->create_3d_T_array(lattice,1,nx_local+2,ny_local+2,
					"app_grain:lattice");
    lat_2d = lattice[0];

    // init local and ghost spins to zero

    for (i = 0; i <= nx_local+1; i++) 
      for (j = 0; j <= ny_local+1; j++) 
	lat_2d[i][j] = 0;
  
    // initialize local spins
  
    random = new RandomPark(seed);
    for (i = 0; i < 100; i++) random->uniform();
  
    // loop over global list
    // so that assigment is independent of parallel decomposition
    // and also so that each local domain is initialized with
    // different spins.
    for (i = 1; i <= nx_global; i++) {
      ii = i - nx_offset;
      if (ii >= 1 && ii <= nx_local) { 
	for (j = 1; j <= ny_global; j++) {
	  jj = j - ny_offset;
	  if (jj >= 1 && jj <= ny_local) { 
	    lat_2d[ii][jj] = random->irandom(nspins);
	  } else {
	    random->irandom(nspins);
	  }
	}
      } else {
	for (j = 1; j <= ny_global; j++) {
	  random->irandom(nspins);
	}
      }
    }

  } else if (dimension == 3) {

    // partition lattice across procs
    // allocate my sub-section
  
    procs2lattice_3d();
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

    if (nx_local < 2 || ny_local < 2 || nz_local < 2)
      error->one("Lattice per proc is too small");

    // Figure out who neighbor processors are    

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

    memory->create_3d_T_array(lattice,nx_local+2,ny_local+2,nz_local+2,
					"app_grain:lattice");

    lat_3d = lattice;

    // init local and ghost spins to zero

    for (i = 0; i <= nx_local+1; i++) 
      for (j = 0; j <= ny_local+1; j++) 
	for (k = 0; k <= nz_local+1; k++) 
	  lat_3d[i][j][k] = 0;
  
    // initialize local spins
  
    random = new RandomPark(seed);
    for (i = 0; i < 100; i++) random->uniform();
  
    // loop over global list
    // so that assigment is independent of parallel decomposition
    // and also so that each local domain is initialized with
    // different spins.
    for (i = 1; i <= nx_global; i++) {
      ii = i - nx_offset;
      if (ii >= 1 && ii <= nx_local) { 
	for (j = 1; j <= ny_global; j++) {
	  jj = j - ny_offset;
	  if (jj >= 1 && jj <= ny_local) { 
	    for (k = 1; k <= nz_global; k++) {
	      kk = k - nz_offset;
	      if (kk >= 1 && kk <= nz_local) { 
		lat_3d[ii][jj][kk] = random->irandom(nspins);
	      } else {
		random->irandom(nspins);
	      }
	    }
	  } else {
	    for (k = 1; k <= nz_global; k++) {
	      random->irandom(nspins);
	    }
	  }
	}
      } else {
	for (j = 1; j <= ny_global; j++) {
	  for (k = 1; k <= nz_global; k++) {
	    random->irandom(nspins);
	  }
	}
      }
    }
  }
  
  // setup other classes
  
  ((SweepGrain*)sweep)->init(this, me, nprocs, dimension, 
	      nx_local, ny_local, nz_local, 
	      nx_global, ny_global, nz_global,  
	      nx_offset, ny_offset, nz_offset, 
	      lattice,
	      procwest, proceast, 
	      procsouth, procnorth, 
	      procdown, procup, 
	      nspins, temperature);

  // Print layout info
  
  if (me == 0) {
    if (screen) {
      fprintf(screen," nx_global = %d \n ny_global = %d \n nz_global = %d \n",
            nx_global,ny_global,nz_global);
      fprintf(screen," nx_procs = %d \n ny_procs = %d \n nz_procs = %d \n",
            nx_procs,ny_procs,nz_procs);
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

//   char *fstring = "Iteration %6d, in dump()";
//   int len = strlen(fstring)+32;
//   char *title = new char[len];
//   sprintf(title,fstring,0);
//   dump_detailed(title);
//   delete [] title;

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
  if (sweep == NULL) error->all("No sweep class defined");
  
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

void AppGrain::iterate()
{
  timer->barrier_start(TIME_LOOP);
  
  for (int i = 0; i < nsweep; i++) {
    ntimestep++;

    sweep->do_sweep();

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

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppGrain::stats()
{
  double energy;
  
  energy = sweep->compute_energy();
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

  delete [] dumpbuf;
  maxdumpbuf = 0;
  int mybuf = 4*nx_local*ny_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);
  dumpbuf = new int[maxdumpbuf];
  
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
      dumpbuf[m++] = 2*n + 1;
      dumpbuf[m++] = lat_2d[i][j];
      dumpbuf[m++] = 2*n + 2;
      dumpbuf[m++] = lat_2d[i][j];
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

//   // Dump spins in dump
//   char *fstring = "Iteration %6d, in dump()";
//   int len = strlen(fstring)+32;
//   char *title = new char[len];
//   sprintf(title,fstring,ntimestep);
//   dump_detailed(title);
//   delete [] title;

}

/* ----------------------------------------------------------------------
   dump a snapshot of grains to the screen in 2D layout
   all the spins for each processor domain are printed out
------------------------------------------------------------------------- */

void AppGrain::dump_detailed(char* title)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxhtmp,nyhtmp,nzhtmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff
  if (dimension == 2) {
    maxbuftmp = ((nx_global-1)/nx_procs+3)*((ny_global-1)/ny_procs+3)+9;
    nsend = (nx_local+2)*(ny_local+2)+9;
    if (maxbuftmp < nsend) 
      error->one("maxbuftmp size too small in AppGrain::dump_detailed()");
  } else { 
    maxbuftmp = ((nx_global-1)/nx_procs+3)*((ny_global-1)/ny_procs+3)*((nz_global-1)/nz_procs+3)+9;
    nsend = (nx_local+2)*(ny_local+2)*(nz_local+2)+9;
    if (maxbuftmp < nsend)
      error->one("maxbuftmp size too small in AppGrain::dump_detailed()");
  }

  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"appgrain:dump_detailed:buftmp");

  // proc 0 writes interactive dump header

  if (me == 0) {
    if (screen) {
      fprintf(screen,"*** Interactive Dump ***\n");
      fprintf(screen,"Title = %s\n",title);
      fprintf(screen,"nx_global = %d ny_global = %d nz_global = %d \n",nx_global,ny_global,nz_global);
    }
  }

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = nz_local;
  // Need to delete these two
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = nz_offset;

  // pack my lattice values into buffer
  // Need to violate normal ordering in order to simplify output

  if (dimension == 2) {
    for (int j = 0; j <= ny_local+1; j++) {
      for (int i = 0; i <= nx_local+1; i++) {
	buftmp[m++] = lat_2d[i][j];
      }
    }
  } else {
    for (int k = 0; k <= nz_local+1; k++) {
      for (int j = 0; j <= ny_local+1; j++) {
	for (int i = 0; i <= nx_local+1; i++) {
	  buftmp[m++] = lat_3d[i][j][k];
	}
      }
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
	fprintf(screen,"nxlocal = %d \nnylocal = %d \nnzlocal = %d \n",nxtmp,nytmp,nztmp);
	// Need to delete this one
	fprintf(screen,"0 = %d \n0 = %d \n0 = %d \n",nxhtmp,nyhtmp,nzhtmp);
	fprintf(screen,"nx_offset = %d \nny_offset = %d \nnz_offset = %d \n",nxotmp,nyotmp,nzotmp);
	if (dimension == 2) {
	  m = nrecv;
	  for (int j = nytmp+1; j >= 0; j--) {
	    m-=nxtmp+2;
	    for (int i = 0; i <= nxtmp+1; i++) {
	      fprintf(screen,"%3d",buftmp[m++]);
	    }
	    fprintf(screen,"\n");
	    m-=nxtmp+2;
	  }

	} else {
	  m = nrecv;
	  for (int k = nztmp+1; k >= 0; k--) {
	    fprintf(screen,"\n*** k = %d *** \n",k);
	    for (int j = nytmp+1; j >= 0; j--) {
	      m-=nxtmp+2;
	      for (int i = 0; i <= nxtmp+1; i++) {
		fprintf(screen,"%3d",buftmp[m++]);
	      }
	      fprintf(screen,"\n");
	      m-=nxtmp+2;
	    }
	  }
	}
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
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

void AppGrain::procs2lattice_2d()
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
   assign nprocs to 3d nx/ny/nzglobal lattice so as to minimize surface 
------------------------------------------------------------------------- */

void AppGrain::procs2lattice_3d()
{
  int ipx,ipy,ipz,nremain;
  // Need float arithmetic to handle uneven ratios
  double boxx,boxy,boxz,surf;
  double bestsurf = 2 * (nx_global*ny_global + ny_global*nz_global + 
                      nz_global*nx_global);
  
  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // for 2d, insure ipz = 1

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
}

int AppGrain::energy_site(int ik, int i, int j) {
  int nk;

  nk = 0;
  if (ik != lat_2d[i-1][j-1]) nk++;
  if (ik != lat_2d[i-1][j]) nk++;
  if (ik != lat_2d[i-1][j+1]) nk++;
  if (ik != lat_2d[i][j-1]) nk++;
  if (ik != lat_2d[i][j+1]) nk++;
  if (ik != lat_2d[i+1][j-1]) nk++;
  if (ik != lat_2d[i+1][j]) nk++;
  if (ik != lat_2d[i+1][j+1]) nk++;

  return nk;
}

void AppGrain::update_mask(char** mask, int i, int j) {

  // Unset masks for self and neighbors
  mask[i-1][j-1] = 0;
  mask[i-1][j] = 0;
  mask[i-1][j+1] = 0;
  mask[i][j-1] = 0;
  mask[i][j] = 0;
  mask[i][j+1] = 0;
  mask[i+1][j-1] = 0;
  mask[i+1][j] = 0;
  mask[i+1][j+1] = 0;
}

int AppGrain::energy_site(int ik, int i, int j, int k) {
  int nk;

  nk = 0;
  if (ik != lat_3d[i-1][j-1][k-1]) nk++;
  if (ik != lat_3d[i-1][j-1][k]) nk++;
  if (ik != lat_3d[i-1][j-1][k+1]) nk++;

  if (ik != lat_3d[i-1][j][k-1]) nk++;
  if (ik != lat_3d[i-1][j][k]) nk++;
  if (ik != lat_3d[i-1][j][k+1]) nk++;

  if (ik != lat_3d[i-1][j+1][k-1]) nk++;
  if (ik != lat_3d[i-1][j+1][k]) nk++;
  if (ik != lat_3d[i-1][j+1][k+1]) nk++;

  if (ik != lat_3d[i][j-1][k-1]) nk++;
  if (ik != lat_3d[i][j-1][k]) nk++;
  if (ik != lat_3d[i][j-1][k+1]) nk++;

  if (ik != lat_3d[i][j][k-1]) nk++;
  if (ik != lat_3d[i][j][k+1]) nk++;

  if (ik != lat_3d[i][j+1][k-1]) nk++;
  if (ik != lat_3d[i][j+1][k]) nk++;
  if (ik != lat_3d[i][j+1][k+1]) nk++;

  if (ik != lat_3d[i+1][j-1][k-1]) nk++;
  if (ik != lat_3d[i+1][j-1][k]) nk++;
  if (ik != lat_3d[i+1][j-1][k+1]) nk++;

  if (ik != lat_3d[i+1][j][k-1]) nk++;
  if (ik != lat_3d[i+1][j][k]) nk++;
  if (ik != lat_3d[i+1][j][k+1]) nk++;

  if (ik != lat_3d[i+1][j+1][k-1]) nk++;
  if (ik != lat_3d[i+1][j+1][k]) nk++;
  if (ik != lat_3d[i+1][j+1][k+1]) nk++;

  return nk;
}

void AppGrain::update_mask(char*** mask, int i, int j, int k) {

  // Unset masks for self and neighbors
  mask[i-1][j-1][k-1] = 0;
  mask[i-1][j-1][k] = 0;
  mask[i-1][j-1][k+1] = 0;

  mask[i-1][j][k-1] = 0;
  mask[i-1][j][k] = 0;
  mask[i-1][j][k+1] = 0;

  mask[i-1][j+1][k-1] = 0;
  mask[i-1][j+1][k] = 0;
  mask[i-1][j+1][k+1] = 0;

  mask[i][j-1][k-1] = 0;
  mask[i][j-1][k] = 0;
  mask[i][j-1][k+1] = 0;

  mask[i][j][k-1] = 0;
  mask[i][j][k] = 0;
  mask[i][j][k+1] = 0;

  mask[i][j+1][k-1] = 0;
  mask[i][j+1][k] = 0;
  mask[i][j+1][k+1] = 0;

  mask[i+1][j-1][k-1] = 0;
  mask[i+1][j-1][k] = 0;
  mask[i+1][j-1][k+1] = 0;

  mask[i+1][j][k-1] = 0;
  mask[i+1][j][k] = 0;
  mask[i+1][j][k+1] = 0;

  mask[i+1][j+1][k-1] = 0;
  mask[i+1][j+1][k] = 0;
  mask[i+1][j+1][k+1] = 0;
}
