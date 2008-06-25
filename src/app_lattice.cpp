/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic ulator
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
#include "output.h"

#include <map>

using namespace SPPARKS_NS;

enum{NONE,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       RANDOM_2D,RANDOM_3D,FILENAME};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 100
#define MAXLINE 256

/* ---------------------------------------------------------------------- */

AppLattice::AppLattice(SPPARKS *spk, int narg, char **arg) : App(spk,narg,arg)
{
  appclass = LATTICE;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // default settings

  ntimestep = 0;
  time = 0.0;
  temperature = 0.0;
  maxdumpbuf = 0;
  dbuf = NULL;
  fp = NULL;
  propensity = NULL;
  site2i = NULL;
  i2site = NULL;

  lattice = NULL;
  iarray = NULL;
  darray = NULL;
  ninteger = ndouble = 0;
  onesite.ivalue = NULL;
  onesite.dvalue = NULL;
  
  // setup communicator for ghost sites

  comm = new CommLattice(spk);

  // app can override these values in its constructor

  dellocal = 0;
  delghost = 1;
}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  delete [] dbuf;
  memory->sfree(propensity);
  memory->sfree(site2i);
  memory->sfree(i2site);

  memory->sfree(lattice);
  for (int i = 0; i < ninteger; i++) memory->sfree(iarray[i]);
  for (int i = 0; i < ndouble; i++) memory->sfree(darray[i]);
  delete [] iarray;
  delete [] darray;

  delete [] onesite.ivalue;
  delete [] onesite.dvalue;

  delete [] latfile;

  memory->sfree(id);
  memory->sfree(owner);
  memory->sfree(index);
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
  sitecustom = 0;

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

    } else if (strcmp(arg[iarg],"site") == 0) {
      if (iarg+3 > narg) error->all("Illegal app_style command");
      ninteger = atoi(arg[iarg+1]);
      ndouble = atoi(arg[iarg+2]);
      if (ninteger == 0 && ndouble == 0) sitecustom = 0;
      else sitecustom = 1;
      if (ninteger) onesite.ivalue = new int[ninteger];
      if (ndouble) onesite.dvalue = new double[ndouble];
      iarg += 3;

    } else error->all("Illegal app_style command ");
  }

  if (latstyle == NONE) error->all("Illegal app_style command");
}

/* ----------------------------------------------------------------------
   generate lattice and processor decomposition
   allocate per-site memory: single lattice value or customized arrays
 ------------------------------------------------------------------------- */

void AppLattice::create_lattice()
{
  if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI || 
      latstyle == SC_6N || latstyle == SC_26N || 
      latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
    structured_lattice();
    ghosts_from_connectivity();
  } else if (latstyle == RANDOM_2D || latstyle == RANDOM_3D) {
    random_lattice();
    ghosts_within_cutoff();
  } else if (latstyle == FILENAME) {
    file_lattice();
    ghosts_from_connectivity();
  }

  if (sitecustom == 0)
    lattice =
      (int *) memory->smalloc((nlocal+nghost)*sizeof(int),"app:lattice");
  else {
    if (ninteger) iarray = new int*[ninteger];
    for (int i = 0; i < ninteger; i++)
      iarray[i] = (int *)
	memory->smalloc((nlocal+nghost)*sizeof(int),"app:iarray");
    if (ndouble) darray = new double*[ndouble];
    for (int i = 0; i < ndouble; i++)
      darray[i] = (double *)
	memory->smalloc((nlocal+nghost)*sizeof(double),"app:darray");
  }

  // DEBUG: connectivity check on distance

  printf("Nglobal,Nlocal,Nghost = %d %d %d\n",nglobal,nlocal,nghost);
  int sum = 0;
  for (int i = 0; i < nlocal; i++) sum += numneigh[i];
  int all;
  MPI_Allreduce(&sum,&all,1,MPI_INT,MPI_SUM,world);
  printf("Total/Ave neighbor connections %d %g\n",all,(double)all/nglobal);

  /*
  for (int i = 0; i < nlocal; i++) {
    printf("Neighs of site %d = %d\n",id[i],numneigh[i]);
    for (int j = 0; j < numneigh[i]; j++) {
      if (neighbor[i][j] < 0 || neighbor[i][j] >= nlocal+nghost) {
	printf("INDICES: %d %d: %d %d %d\n",i,j,nlocal,nghost,nlocal+nghost);
	error->one("Bad neighbor index");
      }
      double dx = xyz[i][0] - xyz[neighbor[i][j]][0];
      double dy = xyz[i][1] - xyz[neighbor[i][j]][1];
      double dz = xyz[i][2] - xyz[neighbor[i][j]][2];
      if (dx > 0.5*xprd) dx -= xprd;
      if (dx < -0.5*xprd) dx += xprd;
      if (dy > 0.5*yprd) dy -= yprd;
      if (dy < -0.5*yprd) dy += yprd;
      if (dz > 0.5*zprd) dz -= zprd;
      if (dz < -0.5*zprd) dz += zprd;
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      printf("  DIST %d %d: %g\n",id[i],id[neighbor[i][j]],r);
    }
  }
  */
}

/* ----------------------------------------------------------------------
   generate structured lattice
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

  // basis sites of each unit cell depend on lattice

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

  if (latstyle == SQ_4N) maxneigh = 4;
  else if (latstyle == SQ_8N) maxneigh = 8;
  else if (latstyle == TRI) maxneigh = 6;
  else if (latstyle == SC_6N) maxneigh = 6;
  else if (latstyle == SC_26N) maxneigh = 26;
  else if (latstyle == FCC) maxneigh = 12;
  else if (latstyle == BCC) maxneigh = 8;
  else if (latstyle == DIAMOND) maxneigh = 4;

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  memory->create_2d_T_array(neighbor,nlocal,maxneigh,"app:neighbor");

  // create connectivity offsets

  memory->create_3d_T_array(cmap,nbasis,maxneigh,4,"app:cmap");
  offsets();

  // generate global lattice connectivity for each site
  // connect() computes global index of Jth neighbor of global site I
  // global index = 1 to Nglobal

  for (i = 0; i < nlocal; i++) {
    numneigh[i] = maxneigh;
    for (j = 0; j < numneigh[i]; j++) {
      neighbor[i][j] = connect(id[i],j);
      if (neighbor[i][j] <= 0 || neighbor[i][j] > nglobal)
	error->all("Bad connectivity result");
    }
  }

  memory->destroy_3d_T_array(cmap);

  // create and initialize other site arrays

  owner = (int *) memory->smalloc(nlocal*sizeof(int),"app:owner");
  index = (int *) memory->smalloc(nlocal*sizeof(int),"app:index");

  for (i = 0; i < nlocal; i++) {
    owner[i] = me;
    index[i] = i;
  }
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

  // create and initialize other site arrays

  owner = (int *) memory->smalloc(nlocal*sizeof(int),"app:owner");
  index = (int *) memory->smalloc(nlocal*sizeof(int),"app:index");

  for (i = 0; i < nlocal; i++) {
    owner[i] = me;
    index[i] = i;
  }
}

/* ----------------------------------------------------------------------
   read lattice from file
   vertices and global connectivity
 ------------------------------------------------------------------------- */

void AppLattice::file_lattice()
{
  char line[MAXLINE];
  char *eof;

  // open file

  if (me == 0) {
    fp = fopen(latfile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open file %s",latfile);
      error->one(str);
    }
  }

  // read header, skip 2 lines, read global values

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);

    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%d",&dimension);
    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%d",&nglobal);
    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%d",&maxneigh);

    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%lg %lg",&boxxlo,&boxxhi);
    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%lg %lg",&boxylo,&boxyhi);
    if (dimension == 3) {
      eof = fgets(line,MAXLINE,fp);
      sscanf(line,"%lg %lg",&boxzlo,&boxzhi);
    } else {
      boxzlo = -0.5;
      boxzhi = 0.5;
    }
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one("Unexpected end of lattice file");
  }

  MPI_Bcast(&dimension,1,MPI_INT,0,world);
  MPI_Bcast(&nglobal,1,MPI_INT,0,world);
  MPI_Bcast(&maxneigh,1,MPI_INT,0,world);
  MPI_Bcast(&boxxlo,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&boxxhi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&boxylo,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&boxyhi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&boxzlo,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&boxzhi,1,MPI_DOUBLE,0,world);

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  // partition domain

  if (dimension == 2) procs2lattice_2d();
  else if (dimension == 3) procs2lattice_3d();

  // read and broadcast list of global vertices
  // keep ones in my sub-domain
  // NOTE: one at a time for now, should be chunked later

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
  }

  id = NULL;
  xyz = NULL;
  int max = 0;
  nlocal = 0;
  
  int idone;
  double xone[3];
  xone[2] = 0.0;

  for (int i = 0; i < nglobal; i++) {
    if (me == 0) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one("Unexpected end of lattice file");
      if (dimension == 2) sscanf(line,"%d %lg %lg",&idone,&xone[0],&xone[1]);
      else sscanf(line,"%d %lg %lg %lg",&idone,&xone[0],&xone[1],&xone[2]);
    }
    MPI_Bcast(&idone,1,MPI_INT,0,world);
    MPI_Bcast(xone,3,MPI_DOUBLE,0,world);

    if (xone[0] < subxlo || xone[0] >= subxhi || 
	xone[1] < subylo || xone[1] >= subyhi || 
	xone[2] < subzlo || xone[2] >= subzhi) continue;

    if (nlocal == max) {
      max += DELTA;
      id = (int *) memory->srealloc(id,max*sizeof(int),"app:id");
      memory->grow_2d_T_array(xyz,max,3,"app:xyz");
    }
    id[nlocal] = idone;
    xyz[nlocal][0] = xone[0];
    xyz[nlocal][1] = xone[1];
    xyz[nlocal][2] = xone[2];
    nlocal++;
  }

  // error check to see if all vertices stored by some proc

  int ntotal;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
  if (ntotal != nglobal) error->all("Vertices read from file incorrectly");

  // allocate arrays to store lattice connectivity

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  memory->create_2d_T_array(neighbor,nlocal,maxneigh,"app:neighbor");

  // read global connectivities for each vertex and store in neighbor
  // read one line at a time, broadcast it, keep it if its my vertex
  // determine if its mine via hash = map of local vertices

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
  }

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;
  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  int j,m,n;
  int neigh[maxneigh];
  char *word;

  for (int i = 0; i < nglobal; i++) {
    if (me == 0) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one("Unexpected end of lattice file");
      idone = atoi(strtok(line," \t\n\r\f"));
      n = 0;
      while (word = strtok(NULL," \t\n\r\f")) neigh[n++] = atoi(word);
    }
    MPI_Bcast(&idone,1,MPI_INT,0,world);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(neigh,n,MPI_INT,0,world);

    loc = hash.find(idone);
    if (loc == hash.end()) continue;
    m = loc->second;
    numneigh[m] = n;
    for (j = 0; j < n; j++) neighbor[m][j] = neigh[j];
  }

  // create and initialize other site arrays

  owner = (int *) memory->smalloc(nlocal*sizeof(int),"app:owner");
  index = (int *) memory->smalloc(nlocal*sizeof(int),"app:index");

  for (int i = 0; i < nlocal; i++) {
    owner[i] = me;
    index[i] = i;
  }
}

/* ----------------------------------------------------------------------
   global connectivity for each owned site is known
   convert to local indices and add ghost sites
 ------------------------------------------------------------------------- */

void AppLattice::ghosts_from_connectivity()
{
  int i,j,k,m,size;

  // put all owned sites in a map
  // key = global ID, value = local index

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;
  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  // for all neighbors of owned sites:
  // check if I own neighbor or it's already in ghost list
  // if neither, add it to ghost list and to map

  int maxghost = 0;
  Ghost *buf = NULL;
  nghost = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < numneigh[i]; j++) {
      m = neighbor[i][j];
      if (hash.find(m) == hash.end()) {
	if (nghost == maxghost) {
	  maxghost += DELTA;
	  buf = (Ghost *) 
	    memory->srealloc(buf,maxghost*sizeof(Ghost),"app:buf");
	}
	buf[nghost].id = m;
	buf[nghost].proc = -1;
	hash.insert(std::pair<int,int> (m,nlocal+nghost));
	nghost++;
      }
    }
  }

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // maxghost = max ghosts on any proc

  int maxsize;
  MPI_Allreduce(&nghost,&maxsize,1,MPI_INT,MPI_MAX,world);

  buf = (Ghost *) memory->srealloc(buf,maxsize*sizeof(Ghost),"app:buf");
  Ghost *bufcopy = (Ghost *) 
    memory->smalloc(maxsize*sizeof(Ghost),"app:bufcopy");

  // cycle ghost list around ring of procs back to self
  // when receive it, fill in info for any sites I own
  // info = me as owning proc, my local index, xyz coords

  MPI_Request request;
  MPI_Status status;

  size = nghost;

  for (int loop = 0; loop < nprocs; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxsize*sizeof(Ghost),MPI_CHAR,prev,0,world,&request);
      MPI_Send(buf,size*sizeof(Ghost),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Ghost);
      memcpy(buf,bufcopy,size*sizeof(Ghost));
    }
    for (i = 0; i < size; i++) {
      if (buf[i].proc >= 0) continue;
      loc = hash.find(buf[i].id);
      if (loc != hash.end() && loc->second < nlocal) {
	buf[i].proc = me;
	buf[i].index = loc->second;
	buf[i].x = xyz[loc->second][0];
	buf[i].y = xyz[loc->second][1];
	buf[i].z = xyz[loc->second][2];
      }
    }
  }

  // reallocate site arrays so can append ghost info

  id = (int *) memory->srealloc(id,(nlocal+nghost)*sizeof(int),"app:id");
  owner = (int *) memory->srealloc(owner,(nlocal+nghost)*sizeof(int),
				   "app:owner");
  index = (int *) memory->srealloc(index,(nlocal+nghost)*sizeof(int),
				   "app:index");
  memory->grow_2d_T_array(xyz,nlocal+nghost,3,"app:xyz");
  
  // original ghost list came back to me around ring
  // extract info for my ghost sites
  // error if any site is not filled in

  for (i = 0; i < nghost; i++) {
    if (buf[i].proc == -1) error->one("Ghost site was not found");
    j = nlocal + i;
    id[j] = buf[i].id;
    owner[j] = buf[i].proc;
    index[j] = buf[i].index;
    xyz[j][0] = buf[i].x;
    xyz[j][1] = buf[i].y;
    xyz[j][2] = buf[i].z;
  }

  // convert all my neighbor connections to local indices

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < numneigh[i]; j++) {
      m = neighbor[i][j];
      loc = hash.find(m);
      if (loc == hash.end()) error->one("Ghost connection was not found");
      neighbor[i][j] = loc->second;
    }

  // clean up

  memory->sfree(buf);
  memory->sfree(bufcopy);
}

/* ----------------------------------------------------------------------
   global connectivity for each owned site is not known
   inferred from cutoff distance
   need to find neighbors, create local connectivity and identify ghost sites
 ------------------------------------------------------------------------- */

void AppLattice::ghosts_within_cutoff()
{
  int i,j;

  // put all owned sites within cutoff of sub-box face into buf

  int maxbuf = 0;
  Ghost *bufsend = NULL;
  int nsend = 0;

  for (i = 0; i < nlocal; i++) {
    if (xyz[i][0] - subxlo <= cutoff || subxhi - xyz[i][0] <= cutoff ||
	xyz[i][1] - subylo <= cutoff || subyhi - xyz[i][1] <= cutoff ||
	xyz[i][2] - subzlo <= cutoff || subzhi - xyz[i][2] <= cutoff) {
      if (nsend == maxbuf) {
	maxbuf += DELTA;
	bufsend = (Ghost *) 
	  memory->srealloc(bufsend,maxbuf*sizeof(Ghost),"app:bufsend");
      }
      bufsend[nsend].id = id[i];
      bufsend[nsend].proc = me;
      bufsend[nsend].index = i;
      bufsend[nsend].x = xyz[i][0];
      bufsend[nsend].y = xyz[i][1];
      bufsend[nsend].z = xyz[i][2];
      nsend++;
    }
  }

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // maxsend = max send sites on any proc

  int maxsize;
  MPI_Allreduce(&nsend,&maxsize,1,MPI_INT,MPI_MAX,world);

  bufsend = (Ghost *) 
    memory->srealloc(bufsend,maxsize*sizeof(Ghost),"app:bufsend");
  Ghost *bufcopy = (Ghost *) 
    memory->smalloc(maxsize*sizeof(Ghost),"app:bufcopy");

  // cycle send list around ring of procs back to self
  // when receive it, extract any sites within cutoff of my sub-box
  // test for within cutoff:
  //   for each dim:
  //     test if site coord or 2 periodic images are between cutoff bounds
  //     all 3 dims must satisfy this criterion to keep site as potential ghost
  // loop < nprocs-1 skips owned sites since never recv them

  MPI_Request request;
  MPI_Status status;

  maxbuf = 0;
  Ghost *bufrecv = NULL;
  int nrecv = 0;

  int flag;
  double coord,coordlo,coordhi;
  int size = nsend;

  for (int loop = 0; loop < nprocs-1; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxsize*sizeof(Ghost),MPI_CHAR,prev,0,world,&request);
      MPI_Send(bufsend,size*sizeof(Ghost),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Ghost);
      memcpy(bufsend,bufcopy,size*sizeof(Ghost));
    }
    for (i = 0; i < size; i++) {
      coord = bufsend[i].x;
      coordlo = bufsend[i].x - xprd;
      coordhi = bufsend[i].x + xprd;
      flag = 0;
      if (coord >= subxlo-cutoff && coord <= subxhi+cutoff) flag = 1;
      if (coordlo >= subxlo-cutoff && coordlo <= subxhi+cutoff) flag = 1;
      if (coordhi >= subxlo-cutoff && coordhi <= subxhi+cutoff) flag = 1;
      if (flag == 0) continue;

      coord = bufsend[i].y;
      coordlo = bufsend[i].y - yprd;
      coordhi = bufsend[i].y + yprd;
      flag = 0;
      if (coord >= subylo-cutoff && coord <= subyhi+cutoff) flag = 1;
      if (coordlo >= subylo-cutoff && coordlo <= subyhi+cutoff) flag = 1;
      if (coordhi >= subylo-cutoff && coordhi <= subyhi+cutoff) flag = 1;
      if (flag == 0) continue;

      coord = bufsend[i].z;
      coordlo = bufsend[i].z - zprd;
      coordhi = bufsend[i].z + zprd;
      flag = 0;
      if (coord >= subzlo-cutoff && coord <= subzhi+cutoff) flag = 1;
      if (coordlo >= subzlo-cutoff && coordlo <= subzhi+cutoff) flag = 1;
      if (coordhi >= subzlo-cutoff && coordhi <= subzhi+cutoff) flag = 1;
      if (flag == 0) continue;

      if (nrecv == maxbuf) {
	maxbuf += DELTA;
	bufrecv = (Ghost *) memory->srealloc(bufrecv,maxbuf*sizeof(Ghost),
					     "app:bufrecv");
      }
      bufrecv[nrecv].id = bufsend[i].id;
      bufrecv[nrecv].proc = bufsend[i].proc;
      bufrecv[nrecv].index = bufsend[i].index;
      bufrecv[nrecv].x = bufsend[i].x;
      bufrecv[nrecv].y = bufsend[i].y;
      bufrecv[nrecv].z = bufsend[i].z;
      nrecv++;
    }
  }

  // count max neighbors thru expensive N^2 loop
  // NOTE: would be faster to bin my owned + ghost sites
  // loop over owned sites and received sites (possible ghosts)
  // each time a received site is within cutoff, it becomes a ghost site
  // increment ghost count if ghost ID not in hash

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  nghost = 0;
  std::map<int,int> hash;

  double dx,dy,dz,rsq;
  double cutsq = cutoff*cutoff;

  for (i = 0; i < nlocal; i++) {
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

    for (j = 0; j < nrecv; j++) {
      dx = xyz[i][0] - bufrecv[j].x;
      dy = xyz[i][1] - bufrecv[j].y;
      dz = xyz[i][2] - bufrecv[j].z;

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
	if (hash.find(bufrecv[j].id) == hash.end()) {
	  hash.insert(std::pair<int,int> (bufrecv[j].id,0));
	  nghost++;
	}
      }
    }
  }

  // reallocate site arrays so can append ghost info

  id = (int *) memory->srealloc(id,(nlocal+nghost)*sizeof(int),"app:id");
  owner = (int *) memory->srealloc(owner,(nlocal+nghost)*sizeof(int),"app:id");
  index = (int *) memory->srealloc(index,(nlocal+nghost)*sizeof(int),"app:id");
  memory->grow_2d_T_array(xyz,nlocal+nghost,3,"app:xyz");

  // allocate arrays to store lattice connectivity

  int tmp = 0;
  for (i = 0; i < nlocal; i++) tmp = MAX(tmp,numneigh[i]);
  MPI_Allreduce(&tmp,&maxneigh,1,MPI_INT,MPI_MAX,world);
  if (maxneigh == 0) error->all("Random lattice has no connectivity");

  memory->create_2d_T_array(neighbor,nlocal,maxneigh,"app:neighbor");

  // generate site connectivity thru expensive N^2 loop
  // NOTE: would be faster to bin my owned + ghost sites
  // loop over owned sites and received sites (possible ghosts)
  // each time a received site is within cutoff, it becomes a ghost site
  // increment ghost count if ghost ID not in hash

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  nghost = 0;
  hash.clear();

  for (i = 0; i < nlocal; i++) {
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

    for (j = 0; j < nrecv; j++) {
      dx = xyz[i][0] - bufrecv[j].x;
      dy = xyz[i][1] - bufrecv[j].y;
      dz = xyz[i][2] - bufrecv[j].z;

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
	if (hash.find(bufrecv[j].id) == hash.end()) {
	  hash.insert(std::pair<int,int> (bufrecv[j].id,0));

	  id[nlocal+nghost] = bufrecv[j].id;
	  owner[nlocal+nghost] = bufrecv[j].proc;
	  index[nlocal+nghost] = bufrecv[j].index;
	  xyz[nlocal+nghost][0] = bufrecv[j].x;
	  xyz[nlocal+nghost][1] = bufrecv[j].y;
	  xyz[nlocal+nghost][2] = bufrecv[j].z;
	  nghost++;
	}
      }
    }
  }

  // clean up

  memory->sfree(bufsend);
  memory->sfree(bufcopy);
  memory->sfree(bufrecv);
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  int i,j,m;

  // app-specific initialization

  init_app();

  // comm init
  
  comm->init(NULL,delghost,dellocal,NULL);

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
  // KMC sweeper will allocate own propensity and site2i arrays

  memory->sfree(propensity);
  memory->sfree(site2i);
  memory->sfree(i2site);

  int nsites = nlocal;

  if (sweep == NULL) {
    propensity = (double*) memory->smalloc(nsites*sizeof(double),
					   "applattice:propensity");
    site2i = (int *) memory->smalloc(nsites*sizeof(int),"applattice:site2i");
  } else {
    propensity = NULL;
    site2i = NULL;
  }

  i2site = (int *) memory->smalloc(nsites*sizeof(int),"applattice:i2site");

  // initialize lattice <-> site mapping arrays
  // they map proc's entire 1d sub-domain to 1d sites and vice versa
  // KMC sweeper will create sector-specific i2site values and ignore app's
  // KMC sweeper will create sector-specific site2i values and ignore app's

  for (i = 0 ; i < nlocal; i++) i2site[i] = i;
  if (site2i) for (i = 0; i < nsites; i++) site2i[i] = i;

  // initialize sweeper
  
  if (sweep) sweep->init();

  // initialize propensities for solver
  // if KMC sweep, sweeper does its own init of its propensity arrays

  if (propensity) {
    comm->all();

    for (i = 0 ; i < nlocal; i++)
      propensity[i] = site_propensity(i,0);

    solve->init(nlocal,propensity);
  }

  // Initialize output

  output->init(time);

}

/* ---------------------------------------------------------------------- */

void AppLattice::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) output->set_dump(narg,arg);
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
  int i,isite;
  double dt;

  if (time >= stoptime) return;

  timer->barrier_start(TIME_LOOP);
  
  int done = 0;
  
  while (!done) {

    if (propensity) {
      timer->stamp();
      isite = solve->event(&dt);
      timer->stamp(TIME_SOLVE);
      
      if (isite < 0) done = 1;
      else {
	ntimestep++;
	i = site2i[isite];
	site_event(i,1);
	time += dt;
	timer->stamp(TIME_APP);
      }

    } else {
      sweep->do_sweep(dt);
      time += dt;
    }

    if (time >= stoptime) done = 1;

    // Do output

    output->compute(time,done);
    timer->stamp(TIME_OUTPUT);

  }
  
  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice::stats(char *strtmp)
{
  int ntimestepall;
  MPI_Allreduce(&ntimestep,&ntimestepall,1,MPI_INT,MPI_SUM,world);
  sprintf(strtmp," %10d %10g",ntimestepall,time);
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s","Step","Time");
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

void AppLattice::set_dump(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal dump command");

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

  subxlo = boxxlo + iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = boxxlo + (iprocx+1) * xprd/nx_procs;
  else subxhi = boxxhi;

  subylo = boxylo + iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = boxylo + (iprocy+1) * yprd/ny_procs;
  else subyhi = boxyhi;

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

  subxlo = boxxlo + iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = boxxlo + (iprocx+1) * xprd/nx_procs;
  else subxhi = boxxhi;

  subylo = boxylo + iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = boxylo + (iprocy+1) * yprd/ny_procs;
  else subyhi = boxyhi;

  subzlo = boxzlo + iprocz * zprd/nz_procs;
  if (iprocz < nz_procs-1) subzhi = boxzlo + (iprocz+1) * zprd/nz_procs;
  else subzhi = boxzhi;
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

    cmap[1][0][0] = -1; cmap[1][0][1] =  0; cmap[1][0][2] = 1;
    cmap[1][1][0] =  1; cmap[1][1][1] =  0; cmap[1][1][2] = 1;
    cmap[1][2][0] =  0; cmap[1][2][1] =  0; cmap[1][2][2] = 0;
    cmap[1][3][0] =  0; cmap[1][3][1] =  1; cmap[1][3][2] = 0;
    cmap[1][4][0] =  1; cmap[1][4][1] =  0; cmap[1][4][2] = 0;
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
    int m = 0;                             // 1st basis site
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

    m = 0;                                 // 2nd basis site
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

    m = 0;                                 // 3rd basis site
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

    m = 0;                                 // 4th basis site
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
      for (m = 0; m < maxneigh; m++) {
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

/* ----------------------------------------------------------------------
   save the state of a site
   only called for sites with general data, not for default lattice
 ------------------------------------------------------------------------- */

void AppLattice::site_save(int i)
{
  int m;
  for (m = 0; m < ninteger; m++) onesite.ivalue[m] = iarray[m][i];
  for (m = 0; m < ndouble; m++) onesite.dvalue[m] = darray[m][i];
}

/* ----------------------------------------------------------------------
   restore the state of a site
   only called for sites with general data, not for default lattice
 ------------------------------------------------------------------------- */

void AppLattice::site_restore(int i)
{
  int m;
  for (m = 0; m < ninteger; m++) iarray[m][i] = onesite.ivalue[m];
  for (m = 0; m < ndouble; m++) darray[m][i] = onesite.dvalue[m];
}

/* ----------------------------------------------------------------------
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */

void AppLattice::push_connected_neighbors(int i, int* cluster_ids, int id, std::stack<int>* cluststack)
{
  int ii;
  int isite = lattice[i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (lattice[ii] == isite && cluster_ids[ii] == 0) {
      cluststack->push(ii);
      cluster_ids[ii] = id;
    }
  }
}

/* ----------------------------------------------------------------------
   Add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */

void AppLattice::connected_ghosts(int i, int* cluster_ids, Cluster* clustlist, int idoffset)
{
  int iclust;
  int ii;
  int isite = lattice[i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (lattice[ii] == isite && ii >= nlocal) {
      iclust = cluster_ids[i]-idoffset;
      // Add ghost cluster to neighbors of local cluster
      clustlist[iclust].add_neigh(cluster_ids[ii]);
    }
  }
}

