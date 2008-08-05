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

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sweep_lattice3d.h"
#include "app_lattice3d.h"
#include "comm_lattice3d.h"
#include "random_park.h"
#include "solve.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define DELTA 100

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SweepLattice3d::SweepLattice3d(SPPARKS *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  
  // parse optional args

  Lmask = false;
  mask = NULL;
  Lstrict = false;
  Lkmc = false;
  Ladapt = false;
  ranlat = NULL;
  delt = 1.0;
  deln0 = 0.0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"strict") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lstrict = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lstrict = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mask") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lmask = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lmask = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kmc") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lkmc = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lkmc = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"adapt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Ladapt = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Ladapt = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      delt = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"deln") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      Ladapt = true;
      deln0 = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal sweep_style command");
  }

  // set sweep method function ptr

  if (Lkmc && (Lmask || Lstrict))
      error->all("Invalid combination of sweep flags");
  if (Ladapt && !Lkmc)
      error->all("Invalid combination of sweep flags");

  if (Lkmc) sweeper = &SweepLattice3d::sweep_sector_kmc;
  else {
    if (!Lmask && !Lstrict) sweeper = &SweepLattice3d::sweep_sector;
    if (Lmask && !Lstrict) sweeper = &SweepLattice3d::sweep_sector_mask;
    if (!Lmask && Lstrict) sweeper = &SweepLattice3d::sweep_sector_strict;
    if (Lmask && Lstrict) sweeper = &SweepLattice3d::sweep_sector_mask_strict;
  }

  // initializations

  nsector = 8;
  for (int i = 0; i < nsector; i++) {
    sector[i].solve = NULL;
    sector[i].site2ijk = NULL;
    sector[i].ijk2site = NULL;
    sector[i].propensity = NULL;
    sector[i].sites = NULL;
    sector[i].border = NULL;
  }

  // communicator needed between sweep sectors

  comm = new CommLattice3d(spk);
}

/* ---------------------------------------------------------------------- */

SweepLattice3d::~SweepLattice3d()
{
  delete random;
  delete comm;

  memory->destroy_3d_T_array(mask,nxlo,nylo,nzlo);
  memory->destroy_3d_T_array(ranlat,1,1,1);

  for (int isector = 0; isector < nsector; isector++) {
    delete sector[isector].solve;
    memory->destroy_2d_T_array(sector[isector].site2ijk);
    memory->destroy_3d_T_array(sector[isector].ijk2site);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->destroy_2d_T_array(sector[isector].border);
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice3d::init()
{
  int i,j,k,m;

  applattice = (AppLattice3d *) app;

  // error check

  if (Lmask && applattice->temperature > 0.0)
    error->all("Cannot mask sweeping with non-zero temperature");

  // get application params

  lattice = applattice->lattice;

  nx_local = applattice->nx_local;
  ny_local = applattice->ny_local;
  nz_local = applattice->nz_local;
  nx_offset = applattice->nx_offset;
  ny_offset = applattice->ny_offset;
  nz_offset = applattice->nz_offset;
  int nx_global = applattice->nx_global;
  int ny_global = applattice->ny_global;
  int nz_global = applattice->nz_global;

  int procwest = applattice->procwest;
  int proceast = applattice->proceast;
  int procsouth = applattice->procsouth;
  int procnorth = applattice->procnorth;
  int procdown = applattice->procdown;
  int procup = applattice->procup;

  // app-specific settings

  delpropensity = applattice->delpropensity;
  delevent = applattice->delevent;
  numrandom = applattice->numrandom;

  nxlo = applattice->nxlo;
  nxhi = applattice->nxhi;
  nylo = applattice->nylo;
  nyhi = applattice->nyhi;
  nzlo = applattice->nzlo;
  nzhi = applattice->nzhi;

  delcol = delpropensity + 1;
  if (Lstrict) ncolor = delcol*delcol*delcol;
  else ncolor = 1;
  
  // setup sectors
  // partition sites by i,j,k values

  for (int isector = 0; isector < nsector; isector++) {
    delete sector[isector].solve;
    memory->destroy_2d_T_array(sector[isector].site2ijk);
    memory->destroy_3d_T_array(sector[isector].ijk2site);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->destroy_2d_T_array(sector[isector].border);
    sector[isector].solve = NULL;
    sector[isector].propensity = NULL;
    sector[isector].site2ijk = NULL;
    sector[isector].ijk2site = NULL;
    sector[isector].sites = NULL;
    sector[isector].border = NULL;
  }

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

  sector[0].xlo = 1;
  sector[0].xhi = nx_half-1;
  sector[0].ylo = 1;
  sector[0].yhi = ny_half-1;
  sector[0].zlo = 1;
  sector[0].zhi = nz_half-1;
  sector[0].nx = sector[0].xhi - sector[0].xlo + 1;
  sector[0].ny = sector[0].yhi - sector[0].ylo + 1;
  sector[0].nz = sector[0].zhi - sector[0].zlo + 1;
  
  sector[1].xlo = 1;
  sector[1].xhi = nx_half-1;
  sector[1].ylo = 1;
  sector[1].yhi = ny_half-1;
  sector[1].zlo = nz_half;
  sector[1].zhi = nz_local;
  sector[1].nx = sector[1].xhi - sector[1].xlo + 1;
  sector[1].ny = sector[1].yhi - sector[1].ylo + 1;
  sector[1].nz = sector[1].zhi - sector[1].zlo + 1;

  sector[2].xlo = 1;
  sector[2].xhi = nx_half-1;
  sector[2].ylo = ny_half;
  sector[2].yhi = ny_local;
  sector[2].zlo = 1;
  sector[2].zhi = nz_half-1;
  sector[2].nx = sector[2].xhi - sector[2].xlo + 1;
  sector[2].ny = sector[2].yhi - sector[2].ylo + 1;
  sector[2].nz = sector[2].zhi - sector[2].zlo + 1;

  sector[3].xlo = 1;
  sector[3].xhi = nx_half-1;
  sector[3].ylo = ny_half;
  sector[3].yhi = ny_local;
  sector[3].zlo = nz_half;
  sector[3].zhi = nz_local;
  sector[3].nx = sector[3].xhi - sector[3].xlo + 1;
  sector[3].ny = sector[3].yhi - sector[3].ylo + 1;
  sector[3].nz = sector[3].zhi - sector[3].zlo + 1;

  sector[4].xlo = nx_half;
  sector[4].xhi = nx_local;
  sector[4].ylo = 1;
  sector[4].yhi = ny_half-1;
  sector[4].zlo = 1;
  sector[4].zhi = nz_half-1;
  sector[4].nx = sector[4].xhi - sector[4].xlo + 1;
  sector[4].ny = sector[4].yhi - sector[4].ylo + 1;
  sector[4].nz = sector[4].zhi - sector[4].zlo + 1;

  sector[5].xlo = nx_half;
  sector[5].xhi = nx_local;
  sector[5].ylo = 1;
  sector[5].yhi = ny_half-1;
  sector[5].zlo = nz_half;
  sector[5].zhi = nz_local;
  sector[5].nx = sector[5].xhi - sector[5].xlo + 1;
  sector[5].ny = sector[5].yhi - sector[5].ylo + 1;
  sector[5].nz = sector[5].zhi - sector[5].zlo + 1;

  sector[6].xlo = nx_half;
  sector[6].xhi = nx_local;
  sector[6].ylo = ny_half;
  sector[6].yhi = ny_local;
  sector[6].zlo = 1;
  sector[6].zhi = nz_half-1;
  sector[6].nx = sector[6].xhi - sector[6].xlo + 1;
  sector[6].ny = sector[6].yhi - sector[6].ylo + 1;
  sector[6].nz = sector[6].zhi - sector[6].zlo + 1;

  sector[7].xlo = nx_half;
  sector[7].xhi = nx_local;
  sector[7].ylo = ny_half;
  sector[7].yhi = ny_local;
  sector[7].zlo = nz_half;
  sector[7].zhi = nz_local;
  sector[7].nx = sector[7].xhi - sector[7].xlo + 1;
  sector[7].ny = sector[7].yhi - sector[7].ylo + 1;
  sector[7].nz = sector[7].zhi - sector[7].zlo + 1;

  // setup site2ijk for each sector

  for (int isector = 0; isector < nsector; isector++) {
    int nsites = sector[isector].nx * sector[isector].ny * sector[isector].nz;
    memory->create_2d_T_array(sector[isector].site2ijk,nsites,3,
			      "sweep3d:site2ij");
    for (m = 0; m < nsites; m++) {
      i = m / sector[isector].ny/sector[isector].nz + 1;
      j = (m / sector[isector].nz) % sector[isector].ny + 1;
      k = m % sector[isector].nz + 1;
      sector[isector].site2ijk[m][0] = i + sector[isector].xlo - 1;
      sector[isector].site2ijk[m][1] = j + sector[isector].ylo - 1;
      sector[isector].site2ijk[m][2] = k + sector[isector].zlo - 1;
    }
  }

  // setup ijk2site for each sector
  // 0 to nsite-1 for owned points in sector, else -1

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++) {
      memory->create_3d_T_array(sector[isector].ijk2site,
				nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
				"sweep3d:ij2site");
      for (i = 1-delpropensity; i <= nx_local+delpropensity; i++) 
	for (j = 1-delpropensity; j <= ny_local+delpropensity; j++)
	  for (k = 1-delpropensity; k <= nz_local+delpropensity; k++)
	    sector[isector].ijk2site[i][j][k] = -1;
      int nsites = sector[isector].nx * sector[isector].ny * 
	sector[isector].nz;
      for (m = 0; m < nsites; m++) {
	i = sector[isector].site2ijk[m][0];
	j = sector[isector].site2ijk[m][1];
	k = sector[isector].site2ijk[m][2];
	sector[isector].ijk2site[i][j][k] = m;
      }
    }
  }

  // setup mask array for owned and ghost values
  // ghost values written to in app::site_event_rejection

  if (Lmask && mask == NULL) {
    memory->create_3d_T_array(mask,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
			      "sweep3d:mask");
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= nz_local; k++)
	  mask[i][j][k] = 0;
  }

  // setup border array used by masking and KMC solver
  // nborder = # of sites in sector influenced by site outside sector
  // border = list of these sites, stored as lattice indices

  if (Lmask || Lkmc) {
    for (int isector = 0; isector < nsector; isector++)
      sector[isector].nborder = find_border_sites(isector);
  }

  // allocate sites array used by KMC solver for border sites

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++)
      sector[isector].sites =
	(int *) memory->smalloc(sector[isector].nborder*sizeof(int),
				"sweep3d:sites");
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_3d_T_array(ranlat,1,nx_local,1,ny_local,1,nz_local,
			      "sweep3d:ranlat");
    int isite;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= nz_local; k++) {
	  isite = (i+nx_offset)*ny_global*nz_global + 
	    (j+ny_offset)*nz_global + k + nz_offset;
	  ranlat[i][j][k].init(seed+isite);
	}
  }

  // init communication for ghost sites

  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup,
	     delpropensity,delevent);

  // setup KMC solver and propensity arrays for each sector
  // propensity init requires all ghost cell info via comm->all()

  if (Lkmc) {
    comm->all(lattice);

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nx * sector[isector].ny * 
	sector[isector].nz;
      sector[isector].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep3d:propensity");
      for (m = 0; m < nsites; m++) {
	i = sector[isector].site2ijk[m][0];
	j = sector[isector].site2ijk[m][1];
	k = sector[isector].site2ijk[m][2];
	sector[isector].propensity[m] = applattice->site_propensity(i,j,k);
      }

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }

  // compute deln0, if not already specified
  // controls future timesteps

  if (Lkmc && Ladapt) {
    pmax = 0.0;
    int ntmp ;
    double ptmp;
    for (int isector = 0; isector < nsector; isector++) {
      ntmp = sector[isector].solve->get_num_active();
      if (ntmp > 0) {
	ptmp = sector[isector].solve->get_total_propensity();
	ptmp /= ntmp;
	pmax = MAX(ptmp,pmax);
      }
    }
    MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    if (deln0 <= 0.0) deln0 = pmaxall*delt;
    else if (pmaxall > 0.0) delt = deln0/pmaxall;
  }
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice3d::do_sweep(double &dt)
{
  if (Ladapt) pmax = 0.0;
  if (!Lkmc) applattice->ntimestep++;

  for (int icolor = 0; icolor < ncolor; icolor++)
    for (int isector = 0; isector < nsector; isector++) {
      timer->stamp();
      comm->sector(lattice,isector);
      timer->stamp(TIME_COMM);

      (this->*sweeper)(icolor,isector);
      timer->stamp(TIME_SOLVE);

      comm->reverse_sector(lattice,isector);
      timer->stamp(TIME_COMM);
    }

  dt = delt;

  // adjust KMC threshold time

  if (Ladapt) {
    MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_SUM,world);
    if (pmaxall > 0.0) delt = deln0/pmaxall;
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector(int icolor, int isector)
{
  int i,j,k;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++)
	applattice->site_event_rejection(i,j,k,random);
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites with masking
   perform events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_mask(int icolor, int isector)
{
  int i,j,k;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  boundary_clear_mask(isector);

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	if (mask[i][j][k]) continue;
	applattice->site_event_rejection(i,j,k,random);
      }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_strict(int icolor, int isector)
{
  int i,j,k,i0,j0,k0;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  i0 = icolor/(delcol*delcol)  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor/delcol % delcol - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  k0 = icolor%delcol  - (nz_offset + zlo-1) % delcol;
  k0 = (k0 < 0) ? k0+delcol : k0;

  for (i = xlo+i0; i <= xhi; i += delcol)
    for (j = ylo+j0; j <= yhi; j += delcol)
      for (k = zlo+k0;  k <= zhi; k += delcol)
	applattice->site_event_rejection(i,j,k,&ranlat[i][j][k]);
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color with masking
   perform events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_mask_strict(int icolor, int isector)
{
  int i,j,k,m,i0,j0,k0;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  boundary_clear_mask(isector);

  i0 = icolor/(delcol*delcol)  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor/delcol % delcol - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  k0 = icolor%delcol  - (nz_offset + zlo-1) % delcol;
  k0 = (k0 < 0) ? k0+delcol : k0;

  for (i = xlo+i0; i <= xhi; i+=delcol)
    for (j = ylo+j0; j <= yhi; j+=delcol)
      for (k = zlo+k0;  k <= zhi; k += delcol) {
	if (mask[i][j][k]) {
	  for (m = 0; m < numrandom; m++) ranlat[i][j][k].uniform();
	  continue;
	}
	applattice->site_event_rejection(i,j,k,&ranlat[i][j][k]);
      }
}

/* ----------------------------------------------------------------------
   generate events in each sector using KMC solver
 ------------------------------------------------------------------------- */

void SweepLattice3d::sweep_sector_kmc(int icolor, int isector)
{
  double dt,time;
  int done,isite,i,j,k;

  // extract sector specific info from octant struct

  int **border = sector[isector].border;
  int nborder = sector[isector].nborder;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int **site2ijk = sector[isector].site2ijk;
  int ***ijk2site = sector[isector].ijk2site;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // solver, propensity, ijk2site array

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  applattice->propensity = propensity;
  applattice->ijk2site = ijk2site;

  // update propensities for sites which neighbor a site outside sector
  // necessary since outside sites may have changed
  // attribute this chunk of time to comm, b/c due to decomposition

  timer->stamp();

  int nsites = 0;

  for (int m = 0; m < nborder; m++) {
    i = border[m][0];
    j = border[m][1];
    k = border[m][2];
    isite = ijk2site[i][j][k];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j,k);
  }

  solve->update(nsites,sites,propensity);
  timer->stamp(TIME_COMM);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    if (isite < 0) done = 1;
    else {
      time += dt;	
      if (time >= delt) done = 1;
      else {
	i = site2ijk[isite][0];
	j = site2ijk[isite][1];
	k = site2ijk[isite][2];
	applattice->site_event(i,j,k,0,random);
	applattice->ntimestep++;
      }
      timer->stamp(TIME_APP);
    }
  }

  // compute maximum sector propensity per site
  // controls future timesteps

  if (Ladapt) {
    int ntmp = applattice->solve->get_num_active();
    if (ntmp > 0) {
      double ptmp = applattice->solve->get_total_propensity();
      ptmp /= ntmp;
      pmax = MAX(ptmp,pmax);
    }
  }
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = NULL;
  applattice->ijk2site = NULL;
}

/* ----------------------------------------------------------------------
   create list of border sites for a sector
   border site = site in sector with a neighbor outside the sector
   neighbor can be another owned site (outside sector) or a ghost
   border = lattice indices of the sites
 ------------------------------------------------------------------------- */

int SweepLattice3d::find_border_sites(int isector)
{
  int **border = NULL;
  int nborder = 0;
  int nmax = 0;

  // loop over all sites in sector
  // only add those within delta of boundary to border

  int i,j,k;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  int delta = delpropensity + delevent;

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	if (i-xlo >= delta && xhi-i >= delta && 
	    j-ylo >= delta && yhi-j >= delta &&
	    k-zlo >= delta && zhi-k >= delta) continue;

	if (nborder == nmax) {
	  nmax += DELTA;
	  memory->grow_2d_T_array(border,nmax,3,"sweep3d:border");
	}
	border[nborder][0] = i;
	border[nborder][1] = j;
	border[nborder][2] = k;
	nborder++;
      }

  sector[isector].border = border;
  return nborder;
}

/* ----------------------------------------------------------------------
   unset all mask values of owned sites in isector whose propensity
     could change due to events on sites one neighbor outside the sector
   border list stores indices of these sites
   their mask value may be out of date, due to state change in other sectors
 ------------------------------------------------------------------------- */

void SweepLattice3d::boundary_clear_mask(int isector)
{
  int **border = sector[isector].border;
  int nborder = sector[isector].nborder;

  for (int m = 0; m < nborder; m++)
    mask[border[m][0]][border[m][1]][border[m][2]] = 0;
}
