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

#include "spktype.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_sites.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "potential.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

// same as in lattice.cpp

enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

enum{BOX,REGION};
enum{DUMMY,IARRAY,DARRAY};

#define DELTALOCAL 10000
#define DELTABUF 10000
#define EPSILON 0.0001

/* ---------------------------------------------------------------------- */

CreateSites::CreateSites(SPPARKS *spk) : Pointers(spk) {}

/* ---------------------------------------------------------------------- */

void CreateSites::command(int narg, char **arg)
{
  if (app == NULL) 
    error->all(FLERR,"Create_sites command before app_style set");
  if (domain->box_exist == 0) 
    error->all(FLERR,"Create_sites command before simulation box is defined");
  if (app->sites_exist == 1) 
    error->all(FLERR,"Cannot create sites after sites already exist");
  if (domain->lattice == NULL)
    error->all(FLERR,"Cannot create sites with undefined lattice");

  if (narg < 1) error->all(FLERR,"Illegal create_sites command");

  int iarg;
  if (strcmp(arg[0],"box") == 0) {
    style = BOX;
    iarg = 1;
  } else if (strcmp(arg[0],"region") == 0) {
    style = REGION;
    if (narg < 2) error->all(FLERR,"Illegal create_sites command");
    nregion = domain->find_region(arg[1]);
    if (nregion == -1) 
      error->all(FLERR,"Create_sites region ID does not exist");
    iarg = 2;
  } else error->all(FLERR,"Illegal create_sites command");

  // parse optional args

  valueflag = DUMMY;
  nbasis = domain->lattice->nbasis;
  basisflag = new int[nbasis+1];
  basis_ivalue = new int[nbasis+1];
  basis_dvalue = new double[nbasis+1];
  for (int i = 1; i <= nbasis; i++) basisflag[i] = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"value") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_sites command");
      valueflag = 1;
      if (strcmp(arg[iarg+1],"site") == 0) {
	valueflag = IARRAY;
	valueindex = 0;
	if (app->iarray == NULL)
	  error->all(FLERR,"Creating a quantity application does not support");
      } else if (arg[iarg+1][0] == 'i') {
	valueflag = IARRAY;
	valueindex = atoi(&arg[iarg+1][1]);
	if (valueindex < 1 || valueindex > app->ninteger)
	  error->all(FLERR,"Creating a quantity application does not support");
	valueindex--;
      } else if (arg[iarg+1][0] == 'd') {
	valueflag = DARRAY;
	valueindex = atoi(&arg[iarg+1][1]);
	if (valueindex < 1 || valueindex > app->ndouble)
	  error->all(FLERR,"Creating a quantity application does not support");
	valueindex--;
      }
      if (valueflag == IARRAY) ivalue = atoi(arg[iarg+2]);
      else dvalue = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_sites command");
      if (valueflag == DUMMY) 
	error->all(FLERR,"Must use value option before basis option "
		   "in create_sites command");
      int ilo,ihi;
      if (nbasis == 0) 
	error->all(FLERR,"Cannot use create_sites basis with random lattice");
      potential->bounds(arg[iarg+1],nbasis,ilo,ihi);
      int count = 0;
      for (int i = ilo; i <= ihi; i++) {
	basisflag[i] = 1;
	if (valueflag == IARRAY) basis_ivalue[i] = atoi(arg[iarg+2]);
	else if (valueflag == DARRAY) basis_dvalue[i] = atof(arg[iarg+2]);
	count++;
      }
      if (count == 0) error->all(FLERR,"Illegal create_sites command");
      iarg += 3;
    } else error->all(FLERR,"Illegal create_sites command");
  }

  // create sites, either on-lattice or off-lattice

  if (domain->me == 0) {
    if (screen) fprintf(screen,"Creating sites ...\n");
    if (logfile) fprintf(logfile,"Creating sites ...\n");
  }

  app->sites_exist = 1;

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  }

  int dimension = domain->dimension;
  latstyle = domain->lattice->style;

  if (latstyle == LINE_2N ||
      latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI || 
      latstyle == SC_6N || latstyle == SC_26N || 
      latstyle == FCC || latstyle == BCC || latstyle == DIAMOND ||
      latstyle == FCC_OCTA_TETRA) {

    xlattice = domain->lattice->xlattice;
    ylattice = domain->lattice->ylattice;
    zlattice = domain->lattice->zlattice;

    structured_lattice();
    if (latticeflag) structured_connectivity();

  } else if (latstyle == RANDOM_1D || latstyle == RANDOM_2D ||
	     latstyle == RANDOM_3D) {
    random_sites();
    if (latticeflag) random_connectivity();
  }

  if (latticeflag) {
    ghosts_from_connectivity(applattice,applattice->delpropensity);
    applattice->print_connectivity();
  }

  delete [] basisflag;
  delete [] basis_ivalue;
  delete [] basis_dvalue;
}

/* ----------------------------------------------------------------------
   generate sites on structured lattice that fits in simulation box
   loop over entire lattice
   if style = REGION, require a site be in region as well
   each proc keeps sites in its sub-domain
 ------------------------------------------------------------------------- */

void CreateSites::structured_lattice()
{
  int dimension = domain->dimension;
  int nonperiodic = domain->nonperiodic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;
  double boxxhi = domain->boxxhi;
  double boxyhi = domain->boxyhi;
  double boxzhi = domain->boxzhi;
  if (dimension <= 1) boxylo = 0.5 * (boxylo+boxyhi);
  if (dimension <= 2) boxzlo = 0.5 * (boxzlo+boxzhi);

  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;

  double **basis = domain->lattice->basis;
  int **iarray = app->iarray;
  double **darray = app->darray;

  // in periodic dims:
  // check that simulation box is integer multiple of lattice spacing
  
  nx = static_cast<int> (domain->xprd / xlattice);
  if (dimension >= 2) ny = static_cast<int> (domain->yprd / ylattice);
  else ny = 1;
  if (dimension == 3) nz = static_cast<int> (domain->zprd / zlattice);
  else nz = 1;

  if (xperiodic && 
      fabs(nx*xlattice - domain->xprd) > EPSILON)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");
  if (dimension > 1 && yperiodic &&
      fabs(ny*ylattice - domain->yprd) > EPSILON)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");
  if (dimension > 2 && zperiodic && 
      fabs(nz*zlattice - domain->zprd) > EPSILON)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");

  // set domain->nx,ny,nz iff style = BOX and system is fully periodic
  // else site IDs may be non-contiguous and/or ordered irregularly

  if (style == BOX && nonperiodic == 0) {
    domain->nx = nx;
    domain->ny = ny;
    domain->nz = nz;
  }

  // if dim is periodic:
  //   lattice origin = lower box boundary
  //   loop bounds = 0 to N-1
  // if dim is non-periodic:
  //   lattice origin = 0.0
  //   loop bounds = enough to tile box completely, with all basis atoms

  if (xperiodic) {
    xorig = boxxlo;
    xlo = 0;
    xhi = nx-1;
  } else {
    xorig = 0.0;
    xlo = static_cast<int> (boxxlo / xlattice);
    while ((xlo+1)*xlattice > boxxlo) xlo--;
    xlo++;
    xhi = static_cast<int> (boxxhi / xlattice);
    while (xhi*xlattice <= boxxhi) xhi++;
    xhi--;
  }

  if (yperiodic) {
    yorig = boxylo;
    ylo = 0;
    yhi = ny-1;
  } else {
    yorig = 0.0;
    ylo = static_cast<int> (boxylo / ylattice);
    while ((ylo+1)*ylattice > boxylo) ylo--;
    ylo++;
    yhi = static_cast<int> (boxyhi / ylattice);
    while (yhi*ylattice <= boxyhi) yhi++;
    yhi--;
  }

  if (zperiodic) {
    zorig = boxzlo;
    zlo = 0;
    zhi = nz-1;
  } else {
    zorig = 0.0;
    zlo = static_cast<int> (boxzlo / zlattice);
    while ((zlo+1)*zlattice > boxzlo) zlo--;
    zlo++;
    zhi = static_cast<int> (boxzhi / zlattice);
    while (zhi*zlattice <= boxzhi) zhi++;
    zhi--;
  }
  
  // generate xyz coords and store them with site ID
  // tile the simulation box from origin, respecting PBC
  // site IDs should be contiguous if style = BOX and fully periodic
  // for non-periodic dims, check if site is within global box
  // for style = REGION, check if site is within region
  // if non-periodic or style = REGION, IDs may not be contiguous

  int i,j,k,m,nlocal;
  double x,y,z;

  int maxlocal = 0;
  siteijk = NULL;

  tagint n = 0;
  for (k = zlo; k <= zhi; k++)
    for (j = ylo; j <= yhi; j++)
      for (i = xlo; i <= xhi; i++)
	for (m = 0; m < nbasis; m++) {
	  n++;
	  x = (i + basis[m][0])*xlattice + xorig;
	  y = (j + basis[m][1])*ylattice + yorig;
	  z = (k + basis[m][2])*zlattice + zorig;

	  if (nonperiodic) {
	    if (!xperiodic && (x < boxxlo || x >= boxxhi)) continue;
	    if (!yperiodic && (y < boxylo || y >= boxyhi)) continue;
	    if (!zperiodic && (z < boxzlo || z >= boxzhi)) continue;
	  }
	  if (style == REGION &&
	      domain->regions[nregion]->match(x,y,z) == 0) continue;

	  if (x < subxlo || x >= subxhi || 
	      y < subylo || y >= subyhi || 
	      z < subzlo || z >= subzhi) continue;

	  if (latticeflag) applattice->add_site(n,x,y,z);
	  else appoff->add_site(n,x,y,z);
	  nlocal = app->nlocal;

	  if (nlocal > maxlocal) {
	    maxlocal += DELTALOCAL;
	    memory->grow(siteijk,maxlocal,4,"create:siteijk");
	  }

	  siteijk[nlocal-1][0] = i;
	  siteijk[nlocal-1][1] = j;
	  siteijk[nlocal-1][2] = k;
	  siteijk[nlocal-1][3] = m;

	  if (valueflag == IARRAY) {
	    if (basisflag[m+1])
	      iarray[valueindex][nlocal-1] = basis_ivalue[m+1];
	    else iarray[valueindex][nlocal-1] = ivalue;
	  } else if (valueflag == DARRAY) {
	    if (basisflag[m+1]) 
	      darray[valueindex][nlocal-1] = basis_dvalue[m+1];
	    else darray[valueindex][nlocal-1] = dvalue;
	  }
	}

  // print site count

  tagint nbig = app->nlocal;
  MPI_Allreduce(&nbig,&app->nglobal,1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (domain->me == 0) {
    if (screen)
      fprintf(screen,"  " TAGINT_FORMAT " sites\n",app->nglobal);
    if (logfile)
      fprintf(logfile,"  " TAGINT_FORMAT " sites\n",app->nglobal);
  }

  // for style = BOX and periodic system, check if nglobal is correct

  if (style == BOX && domain->nonperiodic == 0) {
    nbig = nbasis;
    nbig = nbig*nx*ny*nz;
    if (style == BOX && app->nglobal != nbig)
      error->all(FLERR,"Did not create correct number of sites");
  }
}

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
   respect non-periodic boundaries
   only called for on-lattice models
 ------------------------------------------------------------------------- */

void CreateSites::structured_connectivity()
{
  int i,j,m,max;
  int ineigh,jneigh,kneigh,mneigh;
  tagint gid;
  double xneigh,yneigh,zneigh;

  int nonperiodic = domain->nonperiodic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;
  double boxxhi = domain->boxxhi;
  double boxyhi = domain->boxyhi;
  double boxzhi = domain->boxzhi;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // set maxneigh and allocate idneigh array to store connectivity

  if (latstyle == LINE_2N) maxneigh = 2;
  else if (latstyle == SQ_4N) maxneigh = 4;
  else if (latstyle == SQ_8N) maxneigh = 8;
  else if (latstyle == TRI) maxneigh = 6;
  else if (latstyle == SC_6N) maxneigh = 6;
  else if (latstyle == SC_26N) maxneigh = 26;
  else if (latstyle == FCC) maxneigh = 12;
  else if (latstyle == BCC) maxneigh = 8;
  else if (latstyle == DIAMOND) maxneigh = 4;
  else if (latstyle == FCC_OCTA_TETRA) maxneigh = 26;

  memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

  // create connectivity offsets

  int nbasis = domain->lattice->nbasis;
  double **basis = domain->lattice->basis;
  memory->create(cmap,nbasis,maxneigh,4,"create:cmap");
  offsets(basis);

  // generate global lattice connectivity for each site
  // for non-periodic dims, site must be in global box and not across boundary
  // for style = REGION, check if site is in region
  // FCC_OCTA_TETRA is special case, # of neighs not same for all sites

  tagint nglobal = app->nglobal;
  int nlocal = app->nlocal;
  tagint *id = app->id;
  int *numneigh = applattice->numneigh;

  for (i = 0; i < nlocal; i++) {
    numneigh[i] = 0;
    if (latstyle == FCC_OCTA_TETRA) {
      if ((id[i]-1) % 16 < 8) max = maxneigh;
      else max = 14;
    } else max = maxneigh;
    
    for (j = 0; j < max; j++) {

      // ijkm neigh = indices of neighbor site
      // calculated from siteijk and cmap offsets

      m = siteijk[i][3];
      ineigh = siteijk[i][0] + cmap[m][j][0];
      jneigh = siteijk[i][1] + cmap[m][j][1];
      kneigh = siteijk[i][2] + cmap[m][j][2];
      mneigh = cmap[m][j][3];

      // xyz neigh = coords of neighbor site
      // calculated in same manner that structured_lattice() generated coords

      xneigh = (ineigh + basis[mneigh][0])*xlattice + xorig;
      yneigh = (jneigh + basis[mneigh][1])*ylattice + yorig;
      zneigh = (kneigh + basis[mneigh][2])*zlattice + zorig;

      // remap neighbor coords and indices into periodic box via ijk neigh

      if (xperiodic) {
	if (ineigh < 0) {
	  xneigh += xprd;
	  ineigh += nx;
	}
	if (ineigh >= nx) {
	  xneigh -= xprd;
	  xneigh = MAX(xneigh,boxxlo);
	  ineigh -= nx;
	}
      }
      if (yperiodic) {
	if (jneigh < 0) {
	  yneigh += yprd;
	  jneigh += ny;
	}
	if (jneigh >= ny) {
	  yneigh -= yprd;
	  yneigh = MAX(yneigh,boxylo);
	  jneigh -= ny;
	}
      }
      if (zperiodic) {
	if (kneigh < 0) {
	  zneigh += zprd;
	  kneigh += nz;
	}
	if (kneigh >= nz) {
	  zneigh -= zprd;
	  zneigh = MAX(zneigh,boxzlo);
	  kneigh -= nz;
	}
      }

      // discard neighs that are outside non-periodic box or region

      if (nonperiodic) {
	if (!xperiodic && (xneigh < boxxlo || xneigh >= boxxhi)) continue;
	if (!yperiodic && (yneigh < boxylo || yneigh >= boxyhi)) continue;
	if (!zperiodic && (zneigh < boxzlo || zneigh >= boxzhi)) continue;
      }
      if (style == REGION &&
	  domain->regions[nregion]->match(xneigh,yneigh,zneigh) == 0) continue;

      // gid = global ID of neighbor
      // calculated in same manner that structured_lattice() generated IDs

      gid = (kneigh-zlo)*(yhi-ylo+1)*(xhi-xlo+1)*nbasis + 
	(jneigh-ylo)*(xhi-xlo+1)*nbasis + (ineigh-xlo)*nbasis + mneigh + 1;

      if (style == BOX && nonperiodic == 0 && (gid <= 0 || gid > nglobal))
	error->all(FLERR,"Bad neighbor site ID");

      // add gid to neigh list of site I

      idneigh[i][numneigh[i]++] = gid;
    }
  }

  // delete siteijk and connectivity offsets

  memory->destroy(siteijk);
  memory->destroy(cmap);
}

/* ----------------------------------------------------------------------
   generate random sites
   each proc keeps those in its sub-domain
 ------------------------------------------------------------------------- */

void CreateSites::random_sites()
{
  int dimension = domain->dimension;
  tagint nrandom = domain->lattice->nrandom;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;

  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;

  int **iarray = app->iarray;
  double **darray = app->darray;

  // generate xyz coords and store them with site ID
  // iterate until atom is within region
  // if coords are in my sub-domain, create site

  int nlocal;
  double x,y,z;

  double seed = ranmaster->uniform();
  RandomPark *random = new RandomPark(seed);

  for (tagint n = 1; n <= nrandom; n++) {
    while (1) {
      x = boxxlo + xprd*random->uniform();
      y = boxylo + yprd*random->uniform();
      z = boxzlo + zprd*random->uniform();
      if (dimension < 2) y = 0.0;
      if (dimension < 3) z = 0.0;
      if (style == REGION) {
	if (domain->regions[nregion]->match(x,y,z) == 1) break;
      } else break;
    }

    if (x < subxlo || x >= subxhi || 
	y < subylo || y >= subyhi || 
	z < subzlo || z >= subzhi) continue;
    
    if (latticeflag) applattice->add_site(n,x,y,z);
    else appoff->add_site(n,x,y,z);

    nlocal = app->nlocal;
    if (valueflag == IARRAY) iarray[valueindex][nlocal] = ivalue;
    else if (valueflag == DARRAY) darray[valueindex][nlocal] = dvalue;
  }

  delete random;

  // print site count
  // check if sum of nlocal = nglobal

  tagint nbig = app->nlocal;
  MPI_Allreduce(&nbig,&app->nglobal,1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (domain->me == 0) {
    if (screen)
      fprintf(screen,"  " TAGINT_FORMAT " sites\n",app->nglobal);
    if (logfile)
      fprintf(logfile,"  " TAGINT_FORMAT " sites\n",app->nglobal);
  }

  if (app->nglobal != nrandom)
    error->all(FLERR,"Did not create correct number of sites");
}

/* ----------------------------------------------------------------------
   infer connectivity from neighbors within cutoff distance
   respect non-periodic boundaries
   only called for on-lattice models
 ------------------------------------------------------------------------- */

void CreateSites::random_connectivity()
{
  int i,j;

  int me = domain->me;
  int nprocs = domain->nprocs;
  int dimension = domain->dimension;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  int nlocal = app->nlocal;
  double cutoff = domain->lattice->cutoff;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xhalf = 0.5 * xprd;
  double yhalf = 0.5 * yprd;
  double zhalf = 0.5 * zprd;

  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;

  tagint *id = app->id;
  double **xyz = app->xyz;

  // put all owned sites within cutoff of subdomain face into buf

  int maxbuf = 0;
  Site *bufsend = NULL;
  int nsend = 0;

  for (i = 0; i < nlocal; i++) {
    if (xyz[i][0] - subxlo <= cutoff || subxhi - xyz[i][0] <= cutoff ||
	xyz[i][1] - subylo <= cutoff || subyhi - xyz[i][1] <= cutoff ||
	xyz[i][2] - subzlo <= cutoff || subzhi - xyz[i][2] <= cutoff) {
      if (nsend == maxbuf) {
	maxbuf += DELTABUF;
	bufsend = (Site *) 
	  memory->srealloc(bufsend,maxbuf*sizeof(Site),"create:bufsend");
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
    memory->srealloc(bufsend,maxsize*sizeof(Site),"create:bufsend");
  Site *bufcopy = (Site *) 
    memory->smalloc(maxsize*sizeof(Site),"create:bufcopy");

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
	maxbuf += DELTABUF;
	bufrecv = (Site *) memory->srealloc(bufrecv,maxbuf*sizeof(Site),
					     "create:bufrecv");
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

  int *numneigh = applattice->numneigh;
  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  double dx,dy,dz,rsq;
  double cutsq = cutoff*cutoff;

  for (i = 0; i < nlocal; i++) {
    for (j = i+1; j < nlocal; j++) {
      dx = xyz[i][0] - xyz[j][0];
      dy = xyz[i][1] - xyz[j][1];
      dz = xyz[i][2] - xyz[j][2];

      if (xperiodic && fabs(dx) > xhalf) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
      if (yperiodic && fabs(dy) > yhalf) {
	if (dy < 0.0) dy += yprd;
	else dy -= yprd;
      }
      if (zperiodic && fabs(dz) > zhalf) {
	if (dz < 0.0) dz += zprd;
	else dz -= zprd;
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

      if (xperiodic && fabs(dx) > xhalf) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
      if (yperiodic && fabs(dy) > yhalf) {
	if (dy < 0.0) dy += yprd;
	else dy -= yprd;
      }
      if (zperiodic && fabs(dz) > zhalf) {
	if (dz < 0.0) dz += zprd;
	else dz -= zprd;
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutsq) numneigh[i]++;
    }
  }

  // set maxneigh and allocate idneigh array to store connectivity

  int tmp = 0;
  for (i = 0; i < nlocal; i++) tmp = MAX(tmp,numneigh[i]);
  MPI_Allreduce(&tmp,&maxneigh,1,MPI_INT,MPI_MAX,world);
  if (maxneigh == 0) error->all(FLERR,"Random lattice has no connectivity");

  memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

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

      if (xperiodic && fabs(dx) > xhalf) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
      if (yperiodic && fabs(dy) > yhalf) {
	if (dy < 0.0) dy += yprd;
	else dy -= yprd;
      }
      if (zperiodic && fabs(dz) > zhalf) {
	if (dz < 0.0) dz += zprd;
	else dz -= zprd;
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutsq) {
	idneigh[i][numneigh[i]++] = id[j];
	idneigh[j][numneigh[j]++] = id[i];
      }
    }

    for (j = 0; j < nrecv; j++) {
      dx = xyz[i][0] - bufrecv[j].x;
      dy = xyz[i][1] - bufrecv[j].y;
      dz = xyz[i][2] - bufrecv[j].z;

      if (xperiodic && fabs(dx) > xhalf) {
	if (dx < 0.0) dx += xprd;
	else dx -= xprd;
      }
      if (yperiodic && fabs(dy) > yhalf) {
	if (dy < 0.0) dy += yprd;
	else dy -= yprd;
      }
      if (zperiodic && fabs(dz) > zhalf) {
	if (dz < 0.0) dz += zprd;
	else dz -= zprd;
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq < cutsq) idneigh[i][numneigh[i]++] = bufrecv[j].id;
    }
  }

  // clean up

  memory->sfree(bufsend);
  memory->sfree(bufcopy);
  memory->sfree(bufrecv);
}

/* ----------------------------------------------------------------------
   set maxneigh and initialize idneigh when lattice created via read_sites
   called from read_sites when it reads in sites and neighbors
 ------------------------------------------------------------------------- */

void CreateSites::read_sites(AppLattice *apl)
{
  int i,j;

  maxneigh = apl->maxneigh;
  memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

  int *numneigh = apl->numneigh;
  int **neighbor = apl->neighbor;
  int nlocal = app->nlocal;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < numneigh[i]; j++)
      idneigh[i][j] = neighbor[i][j];
}

/* ----------------------------------------------------------------------
   create ghosts sites around local sub-domain
   only called for on-lattice models
   pass applattice as pointer so can call from ReadSites
   numneigh and global neighbor IDs of each owned site are known as input
   add ghost sites for delpropensity layers
   form neigh list for each layer of ghost sites one layer at a time
   when done, delpropensity-1 layers have a full numneigh and neigh list
     last delpropensity layers has a partial numneigh and neigh list
   convert neighbor IDs from global indices to local indices
 ------------------------------------------------------------------------- */

void CreateSites::ghosts_from_connectivity(AppLattice *apl, int delpropensity)
{
  int i,j,k,m,proc,owner_ghost,index_ghost;
  tagint idglobal,idghost,idrecv;
  double x,y,z;
  tagint *id;
  int *numneigh,**neighbor;
  double **xyz;
  std::map<tagint,int>::iterator loc;
  std::map<tagint,int> hash;

  int me = domain->me;
  int nprocs = domain->nprocs;
  int nlocal = app->nlocal;

  // nchunk = size of one site datum circulated in message

  int nchunk = 7 + maxneigh;

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // loop over delpropensity layers to build up layers of ghosts

  int npreviousghost;
  int nghost = 0;

  for (int ilayer = 0; ilayer < delpropensity; ilayer++) {

    // put all sites (owned + current ghosts) in hash
    // key = global ID, value = local index

    id = app->id;

    hash.clear();
    for (i = 0; i < nlocal+nghost; i++)
      hash.insert(std::pair<tagint,int> (id[i],i));

    // make a list of sites I need
    // loop over neighbors of owned + current ghost sites
    // check if site is already an owned or ghost site or already in list
    // if not, add it to new site list and to hash

    double *buf = NULL;
    int nbuf = 0;
    int maxbuf = 0;
    int nsite = 0;

    numneigh = apl->numneigh;

    for (i = 0; i < nlocal+nghost; i++) {
      for (j = 0; j < numneigh[i]; j++) {
	idglobal = idneigh[i][j];
	if (hash.find(idglobal) == hash.end()) {
	  if (nbuf + nchunk >= maxbuf) {
	    maxbuf += DELTABUF;
	    memory->grow(buf,maxbuf,"create:buf");
	  }
	  buf[nbuf] = idglobal;
	  buf[nbuf+1] = -1;
	  nbuf += nchunk;
	  hash.insert(std::pair<tagint,int> (idglobal,nlocal+nghost+nsite));
	  nsite++;
	}
      }
    }

    // maxsize = max buf size on any proc

    int maxsize;
    MPI_Allreduce(&nbuf,&maxsize,1,MPI_INT,MPI_MAX,world);

    memory->grow(buf,maxsize,"create:buf");
    double *bufcopy;
    memory->create(bufcopy,maxsize,"create:bufcopy");

    // cycle site list around ring of procs back to self
    // when receive it, fill in info for any sites I own
    // info = proc, local index, xyz, numneigh, list of global neighbor IDs
    
    MPI_Request request;
    MPI_Status status;

    xyz = app->xyz;

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
	idrecv = static_cast<tagint> (buf[m++]);
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
	  buf[m++] = idneigh[j][k];
      }
    }

    // original site list came back to me around ring
    // realloc idneigh to store neighbor info for these ghost sites
    // extract info for my new layer of ghost sites
    // reset numneigh after each call to add_ghost() in case realloc occurred
    // error if any site is not filled in

    npreviousghost = nghost;
    nghost += nsite;
    memory->grow(idneigh,nlocal+nghost,maxneigh,"create:idneigh");

    for (i = 0; i < nsite; i++) {
      m = i * nchunk;
      idghost = static_cast<tagint> (buf[m++]);
      owner_ghost = static_cast<int> (buf[m++]);
      if (owner_ghost < 0) error->one(FLERR,"Ghost site was not found");
      index_ghost = static_cast<int> (buf[m++]);
      x = buf[m++];
      y = buf[m++];
      z = buf[m++];

      apl->add_ghost(idghost,x,y,z,owner_ghost,index_ghost);
      numneigh = apl->numneigh;

      j = nlocal + npreviousghost + i;
      numneigh[j] = static_cast<int> (buf[m++]);
      for (k = 0; k < numneigh[j]; k++)
	idneigh[j][k] = static_cast<tagint> (buf[m++]);
    }

    // clean up

    memory->destroy(buf);
    memory->destroy(bufcopy);
  }

  // can now set AppLattice::maxneigh and allocate AppLattice::neighbor

  apl->maxneigh = maxneigh;
  apl->grow(apl->nmax);
  numneigh = apl->numneigh;
  neighbor = apl->neighbor;

  // convert all global neighbors to local indices in AppLattice::neighbor
  // if i is owned or in delpropensity-1 layers, then error if neigh not found
  // if i is ghost in last delpropensity layer, then delete neigh if not found

  for (i = 0; i < nlocal+nghost; i++) {
    j = 0;
    while (j < numneigh[i]) {
      idglobal = idneigh[i][j];
      loc = hash.find(idglobal);
      if (loc != hash.end()) {
	neighbor[i][j] = loc->second;
	j++;
      } else if (i >= nlocal+npreviousghost) {
	numneigh[i]--;
	for (k = j; k < numneigh[i]; k++) idneigh[i][k] = idneigh[i][k+1];
      } else error->one(FLERR,"Ghost connection was not found");
    }
  }

  // no longer need idneigh since AppLattice::neighbor now exists

  memory->destroy(idneigh);
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets(double **basis)
{
  if (latstyle == LINE_2N) {
    cmap[0][0][0] = -1; cmap[0][0][1] = 0; cmap[0][0][2] = 0; cmap[0][0][3] = 0;
    cmap[0][1][0] =  1; cmap[0][1][1] = 0; cmap[0][1][2] = 0; cmap[0][1][3] = 0;
  }

  if (latstyle == SQ_4N)
    for (int m = 0; m < nbasis; m++)
      offsets_2d(m,basis,xlattice,xlattice,maxneigh,cmap[m]);
  else if (latstyle == SQ_8N)
    for (int m = 0; m < nbasis; m++)
      offsets_2d(m,basis,xlattice,sqrt(2.0)*xlattice,maxneigh,cmap[m]);
  else if (latstyle == TRI)
    for (int m = 0; m < nbasis; m++)
      offsets_2d(m,basis,xlattice,xlattice,maxneigh,cmap[m]);

  if (latstyle == SC_6N)
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,xlattice,xlattice,maxneigh,cmap[m]);
  else if (latstyle == SC_26N)
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,xlattice,sqrt(3.0)*xlattice,maxneigh,cmap[m]);
  else if (latstyle == FCC)
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 maxneigh,cmap[m]);
  else if (latstyle == BCC)
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(3.0)/2.0*xlattice,sqrt(3.0)/2.0*xlattice,
		 maxneigh,cmap[m]);
  else if (latstyle == DIAMOND)
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 maxneigh,cmap[m]);

  else if (latstyle == FCC_OCTA_TETRA) {
    for (int m = 0; m < 4; m++) {
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 12,&cmap[m][0]);
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][12]);
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][18]);
    }
    for (int m = 4; m < 8; m++) {
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][0]);
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 12,&cmap[m][6]);
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][18]);
    }
    for (int m = 8; m < nbasis; m++) {
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][0]);
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][8]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets_2d(int ibasis, double **basis, 
			    double cutlo, double cuthi,
			    int ntarget, int **cmapone)
{
  int i,j,m,n;
  double x0,y0,delx,dely,r;

  n = 0;
  x0 = basis[ibasis][0] * xlattice;
  y0 = basis[ibasis][1] * ylattice;
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (m = 0; m < nbasis; m++) {
	delx = (i+basis[m][0])*xlattice - x0;
	dely = (j+basis[m][1])*ylattice - y0;
	r = sqrt(delx*delx + dely*dely);
	if (r > cutlo-EPSILON && r < cuthi+EPSILON) {
	  if (n == ntarget) error->all(FLERR,"Incorrect lattice neighbor count");
	  cmapone[n][0] = i;
	  cmapone[n][1] = j;
	  cmapone[n][2] = 0;
	  cmapone[n][3] = m;
	  n++;
	}
      }
    }
  }

  if (n != ntarget) error->all(FLERR,"Incorrect lattice neighbor count");
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets_3d(int ibasis, double **basis, 
			    double cutlo, double cuthi, 
			    int ntarget, int **cmapone)
{
  int i,j,k,m,n;
  double x0,y0,z0,delx,dely,delz,r;

  n = 0;
  x0 = basis[ibasis][0] * xlattice;
  y0 = basis[ibasis][1] * ylattice;
  z0 = basis[ibasis][2] * zlattice;
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
	for (m = 0; m < nbasis; m++) {
	  delx = (i+basis[m][0])*xlattice - x0;
	  dely = (j+basis[m][1])*ylattice - y0;
	  delz = (k+basis[m][2])*zlattice - z0;
	  r = sqrt(delx*delx + dely*dely + delz*delz);
	  if (r > cutlo-EPSILON && r < cuthi+EPSILON) {
	    if (n == ntarget) error->all(FLERR,"Incorrect lattice neighbor count");
	    cmapone[n][0] = i;
	    cmapone[n][1] = j;
	    cmapone[n][2] = k;
	    cmapone[n][3] = m;
	    n++;
	  }
	}
      }
    }
  }

  if (n != ntarget) error->all(FLERR,"Incorrect lattice neighbor count");
}
