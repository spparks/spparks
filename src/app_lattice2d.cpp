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
  maxdumpbuf = 0;
  ibufread = NULL;
  ibufdump = NULL;
  dbufdump = NULL;
  fp = NULL;
  fpdump = NULL;
  opendxroot = NULL;
  dump_style = COORD;
  propensity = NULL;
  site2ij = NULL;
  ij2site = NULL;

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
  delete [] ibufdump;
  delete [] dbufdump;
  memory->sfree(propensity);
  memory->destroy_2d_T_array(site2ij);
  memory->destroy_2d_T_array(ij2site);

  if (fpdump) {
    fclose(fpdump);
    fpdump = NULL;
  }
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

/* ----------------------------------------------------------------------
   proc 0 writes element connectivity
   one-time info for viz purposes
------------------------------------------------------------------------- */

void AppLattice2d::dump_header()
{
  // setup comm buf for dumping snapshots

  maxdumpbuf = 0;
  int mybuf = nx_local*ny_local;
  MPI_Allreduce(&mybuf,&maxdumpbuf,1,MPI_INT,MPI_MAX,world);

  if (dump_style == COORD) {
    delete [] dbufdump;
    dbufdump = new double[5*maxdumpbuf];
  } else if (dump_style == OPENDX) {
    delete [] ibufdump;
    ibufdump = new int[2*maxdumpbuf];
  } else if (dump_style == LATFILE) {
    delete [] ibufdump;
    ibufdump = new int[2*maxdumpbuf];
  }

  // no header info in file if style = COORD

  if (dump_style == COORD) return;
  if (dump_style == OPENDX) return;

  if (dump_style == LATFILE) {
    int ntimestepall;
    MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

    // proc 0 does one-time write of nodes and element connectivity
    
    if (me) return;
    
    // number nodes: fast in x, slow in y
    
    fprintf(fpdump,"ITEM: TIMESTEP\n");
    fprintf(fpdump,"%d\n",ntimestepall);
    fprintf(fpdump,"ITEM: NUMBER OF NODES\n");
    fprintf(fpdump,"%d\n",(nx_global+1)*(ny_global+1));
    fprintf(fpdump,"ITEM: BOX BOUNDS\n");
    fprintf(fpdump,"%g %g\n",0.0,(double) nx_global);
    fprintf(fpdump,"%g %g\n",0.0,(double) ny_global);
    fprintf(fpdump,"%g %g\n",0.0,0.0);
    fprintf(fpdump,"ITEM: NODES\n");
    
    int i,j;
    int m = 0;
    for (j = 0; j <= ny_global; j++)
      for (i = 0; i <= nx_global; i++) {
	m++;
	fprintf(fpdump,"%d %d %d %d %d\n",m,1,i,j,0);
      }
    
    // number squares: fast in x, slow in y
    // v1,v2,v3,v4 = 4 corner pts of grid cell in counter-clockwise dir
    
    fprintf(fpdump,"ITEM: TIMESTEP\n");
    fprintf(fpdump,"%d\n",ntimestepall);
    fprintf(fpdump,"ITEM: NUMBER OF SQUARES\n");
    fprintf(fpdump,"%d\n",nx_global*ny_global);
    fprintf(fpdump,"ITEM: SQUARES\n");
    
    int v1,v2,v3,v4;
    m = 0;
    for (j = 0; j < ny_global; j++)
      for (i = 0; i < nx_global; i++) {
	v1 = j*(nx_global+1) + i + 1;
	v2 = j*(nx_global+1) + i+1 + 1;
	v3 = (j+1)*(nx_global+1) + i+1 + 1;
	v4 = (j+1)*(nx_global+1) + i + 1;
	m++;
	fprintf(fpdump,"%d %d %d %d %d %d\n",m,1,v1,v2,v3,v4);
      }
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values
   either as lattice sites or as atom coords
------------------------------------------------------------------------- */

void AppLattice2d::dump()
{
  if (dump_style == COORD) {
    dump_coord();
  } else if (dump_style == OPENDX) {
    dump_opendx();
  } else if (dump_style == LATFILE) {
    dump_lattice();
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice coords values, one ATOM per site
   app can provide the coords
------------------------------------------------------------------------- */

void AppLattice2d::dump_coord()
{
  int size_one = 5;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fpdump,"ITEM: TIMESTEP\n");
    fprintf(fpdump,"%d\n",ntimestepall);
    fprintf(fpdump,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fpdump,"%d\n",nx_global*ny_global);

    double boxxlo,boxxhi,boxylo,boxyhi;
    box_bounds(&boxxlo,&boxxhi,&boxylo,&boxyhi);

    fprintf(fpdump,"ITEM: BOX BOUNDS\n");
    fprintf(fpdump,"%g %g\n",boxxlo,boxxhi);
    fprintf(fpdump,"%g %g\n",boxylo,boxyhi);
    fprintf(fpdump,"%g %g\n",0.0,0.0);
    fprintf(fpdump,"ITEM: ATOMS\n");
  }

  // pack my lattice coords into buffer
  // n = global grid cell (0:Nglobal-1)

  double x,y,z;
  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++) {
      n = (ny_offset+j-1)*nx_global + (nx_offset+i-1);
      dbufdump[m++] = n + 1;
      dbufdump[m++] = lattice[i][j];
      xy(i,j,&x,&y);
      dbufdump[m++] = x;
      dbufdump[m++] = y;
      dbufdump[m++] = 0.0;
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
	MPI_Irecv(dbufdump,size_one*maxdumpbuf,MPI_DOUBLE,iproc,0,
		  world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      m = 0;
      for (int i = 0; i < nlines; i++) {
	fprintf(fpdump,"%d %d %g %g %g\n",
		static_cast<int> (dbufdump[m]),static_cast<int> 
		(dbufdump[m+1]),dbufdump[m+2],dbufdump[m+3],dbufdump[m+4]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbufdump,me_size,MPI_DOUBLE,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values, readable by OpenDX
------------------------------------------------------------------------- */

void AppLattice2d::dump_opendx()
{
  int nsites = nx_global*ny_global;
  int* datadx;
  int size_one = 2;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 writes timestep header

  if (me == 0) {
    int lroot = strlen(opendxroot);
    int lnum;
    int lsuf = 3;
    if (opendxcount == 0) {
      lnum = 1;
    } else {
      lnum = int(log(opendxcount)/log(10))+1;
    }
    if (lnum < 5) lnum = 5;
    char filetmp[100];
    if (99 < lroot+lnum+lsuf)
      error->one("app_lattice dump file name too long");
    strcpy(filetmp,opendxroot);
    sprintf(filetmp+lroot,"%05d",opendxcount);
    sprintf(filetmp+lroot+lnum,"%s",".dx");
    fpdump = fopen(filetmp,"w");
    if (!fpdump) error->one("Cannot open app_lattice dump file");
  
    opendxcount++;
    
    fprintf(fpdump,"# app_lattice dump file for OpenDX\n");
    fprintf(fpdump,"# Time = %g\n",time);
    fprintf(fpdump,"# Create regular grid.\n");
    fprintf(fpdump,"object 1 class gridpositions counts %d %d %d\n",
	    nx_global+1,ny_global+1,1+1);
    fprintf(fpdump,"origin  0 0 0 \n");
    fprintf(fpdump,"delta   1 0 0 \n");
    fprintf(fpdump,"delta   0 1 0 \n");
    fprintf(fpdump,"delta   0 0 1 \n");
    fprintf(fpdump,"\n# Create connections.\n");
    fprintf(fpdump,"object 2 class gridconnections counts %d %d %d\n",
	    nx_global+1,ny_global+1,1+1);
    fprintf(fpdump,"\n# Feed data.\n");
    fprintf(fpdump,"object 3 class array type int rank 0 items %d data follows\n#data goes here\n",
	    nsites);
    datadx = (int *) memory->smalloc(nsites*sizeof(int),"diagcluster:datadx");
  }

  // pack my lattice values into buffer
  // n = global grid cell (0:Nglobal-1)

  int n;
  int m = 0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++) {
      n = (ny_offset+j-1)*nx_global + (nx_offset+i-1);
      ibufdump[m++] = n + 1;
      ibufdump[m++] = lattice[i][j];
    }

  int me_size = m;

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;
  int isite;
  
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
	isite = ibufdump[m];
	datadx[isite-1] = ibufdump[m+1];
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(ibufdump,me_size,MPI_INT,0,0,world);
  }

  if (me == 0) {
    isite = 0;
    while (isite < nsites) {
      for (int i = 0; i < 20; i++) {
	fprintf(fpdump,"%d ",datadx[isite++]);
	if (isite == nsites) break;
      }
      fprintf(fpdump,"\n");
    }
    fprintf(fpdump,"attribute \"dep\" string \"connections\"\n");
    fprintf(fpdump,"\n# Create a field.\n");
    fprintf(fpdump,"object \"9 grain microstructure\" class field\n");
    fprintf(fpdump,"component \"positions\" value 1\n");
    fprintf(fpdump,"component \"connections\" value 2\n");
    fprintf(fpdump,"component \"data\" value 3\n");
    fprintf(fpdump,"\nend\n");
    
    fclose(fpdump);
    fpdump = NULL;
    
    memory->sfree(datadx);
  }

}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values, one ELEMENT per site
------------------------------------------------------------------------- */

void AppLattice2d::dump_lattice()
{
  int size_one = 2;

  int ntimestepall;
  MPI_Reduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,0,world);

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fpdump,"ITEM: TIMESTEP\n");
    fprintf(fpdump,"%d\n",ntimestepall);
    fprintf(fpdump,"ITEM: NUMBER OF ELEMENT VALUES\n");
    fprintf(fpdump,"%d\n",nx_global*ny_global);
    fprintf(fpdump,"ITEM: ELEMENT VALUES\n");
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
      ibufdump[m++] = n + 1;
      ibufdump[m++] = lattice[i][j];
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
	fprintf(fpdump,"%d %d\n",ibufdump[m],ibufdump[m+1]);
	m += size_one;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(ibufdump,me_size,MPI_INT,0,0,world);
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
   convert local i,j to x,y
------------------------------------------------------------------------- */

void AppLattice2d::xy(int i, int j, double *x, double *y)
{
  *x = i-1 + nx_offset;
  *y = j-1 + ny_offset;
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

/* ---------------------------------------------------------------------- */

void AppLattice2d::set_dump(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal dump command");

  if (strcmp(arg[1],"coord") == 0) {
    dump_style = COORD;
    if (me == 0) {
      if (fpdump) fclose(fpdump);
      fpdump = fopen(arg[2],"w");
      if (!fpdump) error->one("Cannot open dump file");
    }
  } else if (strcmp(arg[1],"opendx") == 0) {
    dump_style = OPENDX;
    int n = strlen(arg[2]) + 1;
    opendxroot = new char[n];
    strcpy(opendxroot,arg[2]);
    opendxcount = 0;
  } else if (strcmp(arg[1],"lattice") == 0) {
    dump_style = LATFILE;
    if (me == 0) {
      if (fpdump) fclose(fpdump);
      fpdump = fopen(arg[2],"w");
      if (!fpdump) error->one("Cannot open dump file");
    }
  } else error->all("Illegal dump command");
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d global lattice so as to minimize perimeter per proc
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

	ierr = sscanf(line,"%d %d",ipnt,ipnt+1);
	if (ierr != 2) {
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
