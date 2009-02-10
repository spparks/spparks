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
#include "app_lattice.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 100
#define MAXLINE 256
#define CHUNK 1024

/* ---------------------------------------------------------------------- */

void AppLattice::options(int narg, char **arg)
{
  // defaults

  latstyle = NONE;
  latfile = NULL;
  sitecustom = 0;
  infile = NULL;

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

    } else if (strcmp(arg[iarg],"site") == 0) {
      if (iarg+3 > narg) error->all("Illegal app_style command");
      ninteger = atoi(arg[iarg+1]);
      ndouble = atoi(arg[iarg+2]);
      if (ninteger == 0 && ndouble == 0) sitecustom = 0;
      else sitecustom = 1;
      iarg += 3;

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
   generate lattice and processor decomposition
   allocate per-site memory: single lattice value or customized arrays
 ------------------------------------------------------------------------- */

void AppLattice::create_lattice()
{
  if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI || 
      latstyle == SC_6N || latstyle == SC_26N || 
      latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
    structured_lattice();
  } else if (latstyle == RANDOM_2D || latstyle == RANDOM_3D) {
    random_lattice();
    connectivity_within_cutoff();
  } else if (latstyle == FILENAME) {
    file_lattice();
  }

  ghosts_from_connectivity();

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
  int i,n;

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
    if (dimension == 3) z = zprd * random->uniform();
    else z = 0.0;
    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    nlocal++;
  }

  delete random;
  random = new RandomPark(seed);

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

  delete random;

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
  // special treatment to avoid losing vertices at upper boundaries
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

    if (xone[0] < subxlo || xone[1] < subylo || xone[2] < subzlo) continue;
    if (xone[0] > subxhi || xone[1] > subyhi || xone[2] > subzhi) continue;
    if (xone[0] == subxhi && subxhi != boxxhi) continue;
    if (xone[1] == subyhi && subyhi != boxyhi) continue;
    if (xone[2] == subzhi && subzhi != boxzhi) continue;

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

  if (me == 0) fclose(fp);

  // create and initialize other site arrays

  owner = (int *) memory->smalloc(nlocal*sizeof(int),"app:owner");
  index = (int *) memory->smalloc(nlocal*sizeof(int),"app:index");

  for (int i = 0; i < nlocal; i++) {
    owner[i] = me;
    index[i] = i;
  }
}

/* ----------------------------------------------------------------------
   initialize lattice values from a file
 ------------------------------------------------------------------------- */

void AppLattice::read_file()
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
    error->all("Lattice init file has incorrect number of sites");
  if (ninteger == 0 && ndouble == 0) {
    if (nvalues != 1)
      error->all("Lattice init file has incorrect number of values/site");
  } else if (nvalues != ninteger + ndouble)
    error->all("Lattice init file has incorrect number of values/site");

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

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nglobal - nread;
    if (me == 0) {
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of lattice init file");
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
	if (ninteger == 0 && ndouble == 0) lattice[m] = atoi(values[1]);
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

/* ----------------------------------------------------------------------
   inferred connectivity from geometric neighbors within cutoff distance
 ------------------------------------------------------------------------- */

void AppLattice::connectivity_within_cutoff()
{
  int i,j;

  // put all owned sites within cutoff of subdomain face into buf

  int maxbuf = 0;
  Site *bufsend = NULL;
  int nsend = 0;

  for (i = 0; i < nlocal; i++) {
    if (xyz[i][0] - subxlo <= cutoff || subxhi - xyz[i][0] <= cutoff ||
	xyz[i][1] - subylo <= cutoff || subyhi - xyz[i][1] <= cutoff ||
	xyz[i][2] - subzlo <= cutoff || subzhi - xyz[i][2] <= cutoff) {
      if (nsend == maxbuf) {
	maxbuf += DELTA;
	bufsend = (Site *) 
	  memory->srealloc(bufsend,maxbuf*sizeof(Site),"app:bufsend");
      }
      bufsend[nsend].id = id[i];
      bufsend[nsend].proc = me;
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

  bufsend = (Site *) 
    memory->srealloc(bufsend,maxsize*sizeof(Site),"app:bufsend");
  Site *bufcopy = (Site *) 
    memory->smalloc(maxsize*sizeof(Site),"app:bufcopy");

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
  Site *bufrecv = NULL;
  int nrecv = 0;

  int flag;
  double coord,coordlo,coordhi;
  int size = nsend;

  for (int loop = 0; loop < nprocs-1; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxsize*sizeof(Site),MPI_CHAR,prev,0,world,&request);
      MPI_Send(bufsend,size*sizeof(Site),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Site);
      memcpy(bufsend,bufcopy,size*sizeof(Site));
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
	bufrecv = (Site *) memory->srealloc(bufrecv,maxbuf*sizeof(Site),
					     "app:bufrecv");
      }
      bufrecv[nrecv].id = bufsend[i].id;
      bufrecv[nrecv].proc = bufsend[i].proc;
      bufrecv[nrecv].x = bufsend[i].x;
      bufrecv[nrecv].y = bufsend[i].y;
      bufrecv[nrecv].z = bufsend[i].z;
      nrecv++;
    }
  }

  // count max neighbors thru expensive N^2 loop
  // 1st loop over owned sites
  // 2nd loop over owned sites and received sites
  // NOTE: would be faster to bin owned + received sites
  // each time a neighbor is found within cutoff with PBC, increment numneigh

  numneigh = (int *) memory->smalloc(nlocal*sizeof(int),"app:numneigh");
  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

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
      if (rsq < cutsq) numneigh[i]++;
    }
  }

  // allocate neighbor array to store connectivity

  int tmp = 0;
  for (i = 0; i < nlocal; i++) tmp = MAX(tmp,numneigh[i]);
  MPI_Allreduce(&tmp,&maxneigh,1,MPI_INT,MPI_MAX,world);
  if (maxneigh == 0) error->all("Random lattice has no connectivity");

  memory->create_2d_T_array(neighbor,nlocal,maxneigh,"app:neighbor");

  // generate neighbor connectivity thru same expensive N^2 loop
  // 1st loop over owned sites
  // 2nd loop over owned sites and received sites
  // NOTE: would be faster to bin owned + received sites
  // each time a neighbor is found within cutoff with PBC, increment numneigh

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

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
	neighbor[i][numneigh[i]++] = id[j];
	neighbor[j][numneigh[j]++] = id[i];
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
      if (rsq < cutsq) neighbor[i][numneigh[i]++] = bufrecv[j].id;
    }
  }

  // clean up

  memory->sfree(bufsend);
  memory->sfree(bufcopy);
  memory->sfree(bufrecv);
}

/* ----------------------------------------------------------------------
   create ghosts sites around local sub-domain
   numneigh and global neighbor IDs of each owned site are known as input
   add ghost sites for delpropensity layers
   form neigh list for each layer of ghost sites one layer at a time
   when done, delpropensity-1 layers have a full numneigh and neigh list
     last delpropensity layers has a partial numneigh and neigh list
   convert neighbor IDs from global indices to local indices
 ------------------------------------------------------------------------- */

void AppLattice::ghosts_from_connectivity()
{
  int i,j,k,m,idrecv,proc;
  std::map<int,int>::iterator loc;
  std::map<int,int> hash;

  // nchunk = size of one site datum circulated in message

  int nchunk = 7 + maxneigh;

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // loop over delpropensity layers to build up layers of ghosts

  int npreviousghost;
  nghost = 0;

  for (int ilayer = 0; ilayer < delpropensity; ilayer++) {

    // put all sites (owned + current ghosts) in hash
    // key = global ID, value = local index

    hash.clear();
    for (i = 0; i < nlocal+nghost; i++)
      hash.insert(std::pair<int,int> (id[i],i));

    // make a list of sites I need
    // loop over neighbors of owned + current ghost sites
    // check if site is already an owned or ghost site or already in list
    // if not, add it to new site list and to hash

    double *buf = NULL;
    int nbuf = 0;
    int maxbuf = 0;
    int nsite = 0;

    for (i = 0; i < nlocal+nghost; i++) {
      for (j = 0; j < numneigh[i]; j++) {
	m = neighbor[i][j];
	if (hash.find(m) == hash.end()) {
	  if (nbuf + nchunk >= maxbuf) {
	    maxbuf += DELTA;
	    buf = (double *) 
	      memory->srealloc(buf,maxbuf*sizeof(double),"app:buf");
	  }
	  buf[nbuf] = m;
	  buf[nbuf+1] = -1;
	  nbuf += nchunk;
	  hash.insert(std::pair<int,int> (m,nlocal+nghost+nsite));
	  nsite++;
	}
      }
    }

    // maxsize = max buf size on any proc

    int maxsize;
    MPI_Allreduce(&nbuf,&maxsize,1,MPI_INT,MPI_MAX,world);

    buf = (double *) memory->srealloc(buf,maxsize*sizeof(double),"app:buf");
    double *bufcopy = (double *) 
      memory->smalloc(maxsize*sizeof(double),"app:bufcopy");

    // cycle site list around ring of procs back to self
    // when receive it, fill in info for any sites I own
    // info = proc, local index, xyz, numneigh, list of global neighbor IDs
    
    MPI_Request request;
    MPI_Status status;
    
    int size = nbuf;

    for (int loop = 0; loop < nprocs; loop++) {
      if (me != next) {
	MPI_Irecv(bufcopy,maxsize,MPI_DOUBLE,prev,0,world,&request);
	MPI_Send(buf,size,MPI_DOUBLE,next,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&size);
	nsite = size / nchunk;
	memcpy(buf,bufcopy,size*sizeof(double));
      }
      for (int i = 0; i < nsite; i++) {
	m = i * nchunk;
	idrecv = static_cast<int> (buf[m++]);
	proc = static_cast<int> (buf[m++]);
	if (proc >= 0) continue;
	loc = hash.find(idrecv);
	if (loc == hash.end() || loc->second >= nlocal) continue;

	j = loc->second;
	buf[m-1] = me;
	buf[m++] = j;
	buf[m++] = xyz[j][0];
	buf[m++] = xyz[j][1];
	buf[m++] = xyz[j][2];
	buf[m++] = numneigh[j];
	for (k = 0; k < numneigh[j]; k++)
	  buf[m++] = neighbor[j][k];
      }
    }

    // reallocate site arrays so can append next layer of ghosts
    // also reallocate neighbor list arrays

    npreviousghost = nghost;
    nghost += nsite;

    id = (int *) memory->srealloc(id,(nlocal+nghost)*sizeof(int),
				  "app:id");
    owner = (int *) memory->srealloc(owner,(nlocal+nghost)*sizeof(int),
				     "app:owner");
    index = (int *) memory->srealloc(index,(nlocal+nghost)*sizeof(int),
				     "app:index");
    memory->grow_2d_T_array(xyz,nlocal+nghost,3,"app:xyz");

    numneigh = (int *) memory->srealloc(numneigh,(nlocal+nghost)*sizeof(int),
					"app:numneigh");
    memory->grow_2d_T_array(neighbor,nlocal+nghost,maxneigh,"app:neighbor");

    // original site list came back to me around ring
    // extract info for my new layer of ghost sites
    // error if any site is not filled in

    for (i = 0; i < nsite; i++) {
      m = i * nchunk;
      idrecv = static_cast<int> (buf[m++]);
      proc = static_cast<int> (buf[m++]);
      if (proc < 0) error->one("Ghost site was not found");

      j = nlocal + npreviousghost + i;
      id[j] = idrecv;
      owner[j] = proc;
      index[j] = static_cast<int> (buf[m++]);
      xyz[j][0] = buf[m++];
      xyz[j][1] = buf[m++];
      xyz[j][2] = buf[m++];
      numneigh[j] = static_cast<int> (buf[m++]);
      for (k = 0; k < numneigh[j]; k++)
	neighbor[j][k] = static_cast<int> (buf[m++]);
    }

    // clean up

    memory->sfree(buf);
    memory->sfree(bufcopy);
  }

  // convert all my neighbor connections to local indices
  // if i is owned or in delpropensity-1 layers, then error if neigh not found
  // if i is ghost in last delpropensity layer, then delete neigh if not found

  for (i = 0; i < nlocal+nghost; i++) {
    j = 0;
    while (j < numneigh[i]) {
      m = neighbor[i][j];
      loc = hash.find(m);
      if (loc != hash.end()) {
	neighbor[i][j] = loc->second;
	j++;
      } else if (i >= nlocal+npreviousghost) {
	numneigh[i]--;
	for (k = j; k < numneigh[i]; k++) neighbor[i][k] = neighbor[i][k+1];
      } else error->one("Ghost connection was not found");
    }
  }
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
    for (int m = 0; m < maxneigh; m++) cmap[0][m][3] = 0;
 
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
   assign nprocs to global box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void AppLattice::procs2lattice_2d()
{
  int ipx,ipy;
  double boxx,boxy,surf;
  double bestsurf = 2.0 * (xprd+yprd);
  
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
