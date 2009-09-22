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
#include "string.h"
#include "stdlib.h"
#include "app_off_lattice.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

#define DELTA 1000
#define MAXLINE 512
#define CHUNK 1024

/* ---------------------------------------------------------------------- */

void AppOffLattice::options(int narg, char **arg)
{
  // defaults

  latstyle = NONE;
  latfile = NULL;
  px_user = py_user = pz_user = 0;
  px_user = py_user = pz_user = 0;
  infile = NULL;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lattice") == 0) {
      if (iarg+1 > narg) error->all("Illegal app_style command");
      if (strcmp(arg[iarg+1],"line/2n") == 0) latstyle = LINE_2N;
      else if (strcmp(arg[iarg+1],"sq/4n") == 0) latstyle = SQ_4N;
      else if (strcmp(arg[iarg+1],"sq/8n") == 0) latstyle = SQ_8N;
      else if (strcmp(arg[iarg+1],"tri") == 0) latstyle = TRI;
      else if (strcmp(arg[iarg+1],"sc/6n") == 0) latstyle = SC_6N;
      else if (strcmp(arg[iarg+1],"sc/26n") == 0) latstyle = SC_26N;
      else if (strcmp(arg[iarg+1],"fcc") == 0) latstyle = FCC;
      else if (strcmp(arg[iarg+1],"bcc") == 0) latstyle = BCC;
      else if (strcmp(arg[iarg+1],"diamond") == 0) latstyle = DIAMOND;
      else if (strcmp(arg[iarg+1],"fcc/octa/tetra") == 0) 
	latstyle = FCC_OCTA_TETRA;
      else if (strcmp(arg[iarg+1],"random/1d") == 0) latstyle = RANDOM_1D;
      else if (strcmp(arg[iarg+1],"random/2d") == 0) latstyle = RANDOM_2D;
      else if (strcmp(arg[iarg+1],"random/3d") == 0) latstyle = RANDOM_3D;
      else if (strcmp(arg[iarg+1],"file") == 0) latstyle = FILENAME;
      else error->all("Illegal app_style command");

      if (latstyle == LINE_2N) {
	if (iarg+4 > narg) error->all("Illegal app_style command");
	dimension = 1;
	latconst = atof(arg[iarg+2]);
	nx = atoi(arg[iarg+3]);
	iarg += 4;
      } else if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
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
      } else if (latstyle == FCC_OCTA_TETRA) {
	if (iarg+6 > narg) error->all("Illegal app_style command");
	dimension = 3;
	latconst = atof(arg[iarg+2]);
	nx = atoi(arg[iarg+3]);
	ny = atoi(arg[iarg+4]);
	nz = atoi(arg[iarg+5]);
	iarg += 6;
      } else if (latstyle == RANDOM_1D) {
	if (iarg+5 > narg) error->all("Illegal app_style command");
	dimension = 1;
	nrandom = atoi(arg[iarg+2]);
	xprd = atof(arg[iarg+3]);
	cutoff = atof(arg[iarg+4]);
	iarg += 5;
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
	zprd = atof(arg[iarg+5]);
	cutoff = atof(arg[iarg+6]);
	iarg += 7;
      } else if (latstyle == FILENAME) {
	if (iarg+3 > narg) error->all("Illegal app_style command");
	int n = strlen(arg[iarg+2]) + 1;
	latfile = new char[n];
	latfile = strcpy(latfile,arg[iarg+2]);
	iarg += 3;
      }

    } else if (strcmp(arg[iarg],"procs") == 0) {
      if (iarg+4 > narg) error->all("Illegal app_style command");
      px_user = atoi(arg[iarg+1]);
      py_user = atoi(arg[iarg+2]);
      pz_user = atoi(arg[iarg+3]);
      if (px_user || py_user || pz_user) {
	if (px_user*py_user*pz_user != nprocs)
	  error->all("App style proc count is not valid");
      }
      iarg += 4;

    } else if (strcmp(arg[iarg],"input") == 0) {
      if (iarg+2 > narg) error->all("Illegal app_style command");
      int n = strlen(arg[iarg+1]) + 1;
      infile = new char[n];
      strcpy(infile,arg[iarg+1]);
      iarg += 2;

    } else error->all("Illegal app_style command ");
  }

  if (latstyle == NONE) error->all("Illegal app_style command");
}

/* ----------------------------------------------------------------------
   generate processor decomposition and initialize particles
   allocate per-particle memory arrays
 ------------------------------------------------------------------------- */

void AppOffLattice::create_domain()
{
  if (me == 0) {
    if (screen) fprintf(screen,"Creating domain ...\n");
    if (logfile) fprintf(logfile,"Creating domain ...\n");
  }

  if (ninteger) iarray = new int*[ninteger];
  for (int i = 0; i < ninteger; i++) iarray[i] = NULL;
  if (ndouble) darray = new double*[ndouble];
  for (int i = 0; i < ndouble; i++) darray[i] = NULL;

  if (latstyle == LINE_2N ||
      latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI || 
      latstyle == SC_6N || latstyle == SC_26N || 
      latstyle == FCC || latstyle == BCC || latstyle == DIAMOND ||
      latstyle == FCC_OCTA_TETRA) {
    structured_lattice();
  } else if (latstyle == RANDOM_1D || latstyle == RANDOM_2D ||
	     latstyle == RANDOM_3D) {
    random_lattice();
  } else if (latstyle == FILENAME) {
    file_lattice();
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  %d sites\n",nglobal);
    if (logfile) fprintf(logfile,"  %d sites\n",nglobal);
  }
}

/* ----------------------------------------------------------------------
   generate structured lattice
 ------------------------------------------------------------------------- */

void AppOffLattice::structured_lattice()
{
  // determine box extent

  if (latstyle == LINE_2N) {
    boxxlo = 0.0;
    boxxhi = nx * latconst;
    boxylo = -0.5;
    boxyhi = 0.5;
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == SQ_4N || latstyle == SQ_8N) {
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
	     latstyle == FCC || latstyle == BCC || latstyle == DIAMOND ||
	     latstyle == FCC_OCTA_TETRA) {
    boxxlo = boxylo = boxzlo = 0.0;
    boxxhi = nx * latconst;
    boxyhi = ny * latconst;
    boxzhi = nz * latconst;
  }

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  // partition domain

  if (dimension == 1) procs2domain_1d(px_user,py_user,pz_user,
				      xprd,boxxlo,boxxhi,nx_procs,
				      subxlo,subylo,subzlo,
				      subxhi,subyhi,subzhi);
  else if (dimension == 2) procs2domain_2d(px_user,py_user,pz_user,
					   xprd,yprd,
					   boxxlo,boxylo,boxxhi,boxyhi,
					   nx_procs,ny_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);
  else if (dimension == 3) procs2domain_3d(px_user,py_user,pz_user,
					   xprd,yprd,zprd,
					   boxxlo,boxylo,boxzlo,
					   boxxhi,boxyhi,boxzhi,
					   nx_procs,ny_procs,nz_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);

  // basis sites of each unit cell depend on lattice

  int nbasis;
  double **basis;
  memory->create_2d_T_array(basis,16,3,"app:basis");

  if (latstyle == LINE_2N) {
    nbasis = 1;
    nglobal = nbasis * nx;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
  } else if (latstyle == SQ_4N || latstyle == SQ_8N) {
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
  } else if (latstyle == FCC_OCTA_TETRA) {
    nbasis = 16;
    nglobal = nbasis * nx*ny*nz;
    basis[0][0] = 0.0; basis[0][1] = 0.0; basis[0][2] = 0.0;
    basis[1][0] = 0.0; basis[1][1] = 0.5; basis[1][2] = 0.5;
    basis[2][0] = 0.5; basis[2][1] = 0.0; basis[2][2] = 0.5;
    basis[3][0] = 0.5; basis[3][1] = 0.5; basis[3][2] = 0.0;
    basis[4][0] = 0.5; basis[4][1] = 0.0; basis[4][2] = 0.0;
    basis[5][0] = 0.0; basis[5][1] = 0.5; basis[5][2] = 0.0;
    basis[6][0] = 0.0; basis[6][1] = 0.0; basis[6][2] = 0.5;
    basis[7][0] = 0.5; basis[7][1] = 0.5; basis[7][2] = 0.5;
    basis[8][0] = 0.25; basis[8][1] = 0.25; basis[8][2] = 0.25;
    basis[9][0] = 0.75; basis[9][1] = 0.25; basis[9][2] = 0.25;
    basis[10][0] = 0.25; basis[10][1] = 0.75; basis[10][2] = 0.25;
    basis[11][0] = 0.75; basis[11][1] = 0.75; basis[11][2] = 0.25;
    basis[12][0] = 0.25; basis[12][1] = 0.25; basis[12][2] = 0.75;
    basis[13][0] = 0.75; basis[13][1] = 0.25; basis[13][2] = 0.75;
    basis[14][0] = 0.25; basis[14][1] = 0.75; basis[14][2] = 0.75;
    basis[15][0] = 0.75; basis[15][1] = 0.75; basis[15][2] = 0.75;
  }

  if (nglobal <= 0) error->all("Invalid lattice size");

  // generate lattice of sites
  // 1st pass = count lattice points I own in my sub-domain
  // 2nd pass = generate xyz coords and store them with site ID

  if (dimension == 1) ny = nz = 1;
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
	  z = (k + basis[m][2]) * latconstz;

	  if (x < subxlo || x >= subxhi || 
	      y < subylo || y >= subyhi || 
	      z < subzlo || z >= subzhi) continue;

	  nlocal++;
	}

  grow(nlocal);

  nlocal = 0;
  int count = 0;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
	for (m = 0; m < nbasis; m++) {
	  count++;
	  x = (i + basis[m][0]) * latconstx;
	  y = (j + basis[m][1]) * latconsty;
	  z = (k + basis[m][2]) * latconstz;

	  if (x < subxlo || x >= subxhi || 
	      y < subylo || y >= subyhi || 
	      z < subzlo || z >= subzhi) continue;

	  id[nlocal] = count;
	  xyz[nlocal][0] = x;
	  xyz[nlocal][1] = y;
	  xyz[nlocal][2] = z;
	  nlocal++;
	}

  memory->destroy_2d_T_array(basis);
}

/* ----------------------------------------------------------------------
   random lattice
 ------------------------------------------------------------------------- */

void AppOffLattice::random_lattice()
{
  int i,n;

  if (nrandom <= 0) error->all("Invalid lattice size");

  if (latstyle == RANDOM_1D) {
    boxxlo = 0.0;
    boxxhi = xprd;
    boxylo = -0.5;
    boxyhi = 0.5;
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == RANDOM_2D) {
    boxxlo = 0.0;
    boxxhi = xprd;
    boxylo = 0.0;
    boxyhi = yprd;
    boxzlo = -0.5;
    boxzhi = 0.5;
  } else if (latstyle == RANDOM_3D) {
    boxxlo = 0.0;
    boxxhi = xprd;
    boxylo = 0.0;
    boxyhi = yprd;
    boxzlo = 0.0;
    boxzhi = xprd;
  }

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;

  // partition domain

  if (dimension == 1) procs2domain_1d(px_user,py_user,pz_user,
				      xprd,boxxlo,boxxhi,nx_procs,
				      subxlo,subylo,subzlo,
				      subxhi,subyhi,subzhi);
  else if (dimension == 2) procs2domain_2d(px_user,py_user,pz_user,
					   xprd,yprd,
					   boxxlo,boxylo,boxxhi,boxyhi,
					   nx_procs,ny_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);
  else if (dimension == 3) procs2domain_3d(px_user,py_user,pz_user,
					   xprd,yprd,zprd,
					   boxxlo,boxylo,boxzlo,
					   boxxhi,boxyhi,boxzhi,
					   nx_procs,ny_procs,nz_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);

  // generate random sites
  // 1st pass = count sites I own in my sub-domain
  // 2nd pass = generate xyz coords and store them with site ID
  // use new RNG seed twice so as to generate same points

  double seed = ranmaster->uniform();
  RandomPark *random = new RandomPark(seed);

  double x,y,z;
  nglobal = nrandom;
  nlocal = 0;

  for (n = 0; n < nglobal; n++) {
    x = xprd * random->uniform();
    y = yprd * random->uniform();
    z = zprd * random->uniform();
    if (dimension == 1) y = 0.0;
    if (dimension != 3) z = 0.0;
    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    nlocal++;
  }

  delete random;
  random = new RandomPark(seed);

  grow(nlocal);

  nlocal = 0;

  for (n = 0; n < nglobal; n++) {
    x = xprd * random->uniform();
    y = yprd * random->uniform();
    z = zprd * random->uniform();
    if (dimension == 1) y = 0.0;
    if (dimension != 3) z = 0.0;

    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    
    id[nlocal] = n + 1;
    xyz[nlocal][0] = x;
    xyz[nlocal][1] = y;
    xyz[nlocal][2] = z;
    nlocal++;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   read lattice from file
   vertices but no global connectivity
 ------------------------------------------------------------------------- */

void AppOffLattice::file_lattice()
{
  int i,j,m,n,maxneigh;
  FILE *fp;
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

    if (dimension != 1) {
      eof = fgets(line,MAXLINE,fp);
      sscanf(line,"%lg %lg",&boxylo,&boxyhi);
    } else {
      boxylo = -0.5;
      boxyhi = 0.5;
    }

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

  if (dimension == 1) procs2domain_1d(px_user,py_user,pz_user,
				      xprd,boxxlo,boxxhi,nx_procs,
				      subxlo,subylo,subzlo,
				      subxhi,subyhi,subzhi);
  else if (dimension == 2) procs2domain_2d(px_user,py_user,pz_user,
					   xprd,yprd,
					   boxxlo,boxylo,boxxhi,boxyhi,
					   nx_procs,ny_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);
  else if (dimension == 3) procs2domain_3d(px_user,py_user,pz_user,
					   xprd,yprd,zprd,
					   boxxlo,boxylo,boxzlo,
					   boxxhi,boxyhi,boxzhi,
					   nx_procs,ny_procs,nz_procs,
					   subxlo,subylo,subzlo,
					   subxhi,subyhi,subzhi);

  // read and broadcast list of global vertices
  // keep ones in my sub-domain
  // special treatment to avoid losing vertices at upper boundaries
  // include a MPI_Allreduce to force sync, else fast Bcasts crash sometimes

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
  }

  nlocal = 0;
  
  int idone,tmp,tmpall;
  double xone[3];
  xone[1] = xone[2] = 0.0;

  int nchunk;
  int nread = 0;
  char *buffer = new char[CHUNK*MAXLINE];
  char *next,*bufptr;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nglobal - nread;
    if (me == 0) {
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of lattice file");
	if (strlen(&buffer[m]) == MAXLINE-1)
	  error->one("Vertex line too long in lattice file - boost MAXLINE");
	m += strlen(&buffer[m]);
      }
      buffer[m++] = '\n';
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);
    MPI_Allreduce(&tmp,&tmpall,0,MPI_INT,MPI_SUM,world);

    bufptr = buffer;
    for (i = 0; i < nchunk; i++) {
      next = strchr(bufptr,'\n');

      idone = atoi(strtok(bufptr," \t\n\r\f"));
      xone[0] = atof(strtok(NULL," \t\n\r\f"));
      if (dimension != 1) xone[1] = atof(strtok(NULL," \t\n\r\f"));
      if (dimension == 3) xone[2] = atof(strtok(NULL," \t\n\r\f"));

      bufptr = next + 1;

      if (xone[0] < subxlo || xone[1] < subylo || xone[2] < subzlo) continue;
      if (xone[0] > subxhi || xone[1] > subyhi || xone[2] > subzhi) continue;
      if (xone[0] == subxhi && subxhi != boxxhi) continue;
      if (xone[1] == subyhi && subyhi != boxyhi) continue;
      if (xone[2] == subzhi && subzhi != boxzhi) continue;

      if (nlocal == nmax) grow(0);

      id[nlocal] = idone;
      xyz[nlocal][0] = xone[0];
      xyz[nlocal][1] = xone[1];
      xyz[nlocal][2] = xone[2];
      nlocal++;
    }

    nread += nchunk;
  }

  // error check to see if all vertices stored by some proc

  int ntotal;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
  if (ntotal != nglobal) error->all("Vertices read from file incorrectly");

  delete [] buffer;
  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   initialize lattice values from a file
 ------------------------------------------------------------------------- */

void AppOffLattice::read_file()
{
  int i,j,k,m,n;
  int nglobal_file,nvalues;
  FILE *fp;
  char line[MAXLINE];
  char *eof;

  // read header of file

  if (me == 0) {
    fp = fopen(infile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open file %s",infile);
      error->one(str);
    }
  }

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    eof = fgets(line,MAXLINE,fp);
    sscanf(line,"%d %d",&nglobal_file,&nvalues);
    eof = fgets(line,MAXLINE,fp);
  }

  MPI_Bcast(&nglobal_file,1,MPI_INT,0,world);
  MPI_Bcast(&nvalues,1,MPI_INT,0,world);

  if (nglobal_file != nglobal)
    error->all("Input file has incorrect number of sites");
  if (ninteger == 0 && ndouble == 0) {
    if (nvalues != 1)
      error->all("Input file has incorrect number of values/site");
  } else if (nvalues != ninteger + ndouble)
    error->all("Input file has incorrect number of values/site");

  // increment nvalues to include ID

  nvalues++;

  // put all my site IDs into a hash table so can look them up

  std::map<int,int> hash;
  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  // read file in chunks of lines, bcast to other procs
  // tokenize each line into values
  // if I own site ID via hash table lookup, store its values

  int nchunk,idsite;
  int nread = 0;
  char *buffer = new char[CHUNK*MAXLINE];
  char *next,*bufptr;
  char **values = new char*[nvalues];
  std::map<int,int>::iterator loc;

  int site_only;
  if (ninteger == 1 && ndouble == 0) site_only = 1;
  else site_only = 0;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nglobal - nread;
    if (me == 0) {
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of input file");
	m += strlen(&buffer[m]);
      }
      buffer[m++] = '\n';
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    bufptr = buffer;
    for (i = 0; i < nchunk; i++) {
      next = strchr(bufptr,'\n');
      values[0] = strtok(bufptr," \t\n\r\f");
      for (j = 1; j < nvalues; j++)
	values[j] = strtok(NULL," \t\n\r\f");
      idsite = atoi(values[0]);
      loc = hash.find(idsite);
      if (loc != hash.end()) {
	m = loc->second;
	if (site_only) iarray[0][m] = atoi(values[1]);
	else {
	  n = 1;
	  for (k = 0; k < ninteger; k++) iarray[k][m] = atoi(values[n++]);
	  for (k = 0; k < ndouble; k++) darray[k][m] = atof(values[n++]);
	}
      }
      bufptr = next + 1;
    }

    nread += nchunk;
  }

  delete [] values;
  delete [] buffer;

  if (me == 0) fclose(fp);
}
