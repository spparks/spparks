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
#include "sweep_lattice2d.h"
#include "app_lattice2d.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define DELTA 100

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SweepLattice2d::SweepLattice2d(SPPARKS *spk, int narg, char **arg) : 
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
  Ladapt = true;
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
      // This value indicates that delt is specified, not deln
      deln0 = 0.0;
    } else if (strcmp(arg[iarg],"deln") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      deln0 = atof(arg[iarg+1]);
      iarg += 2;
      if (deln0 <= 0.0) error->all("Illegal sweep_style command");
    } else error->all("Illegal sweep_style command");
  }

  // set sweep method function ptr

  if (Lkmc && (Lmask || Lstrict))
      error->all("Invalid combination of sweep flags");

  if (Lkmc) sweeper = &SweepLattice2d::sweep_sector_kmc;
  else {
    if (!Lmask && !Lstrict) sweeper = &SweepLattice2d::sweep_sector;
    if (Lmask && !Lstrict) sweeper = &SweepLattice2d::sweep_sector_mask;
    if (!Lmask && Lstrict) sweeper = &SweepLattice2d::sweep_sector_strict;
    if (Lmask && Lstrict) sweeper = &SweepLattice2d::sweep_sector_mask_strict;
  }

  // initializations

  nsector = 4;
  for (int i = 0; i < nsector; i++) {
    sector[i].solve = NULL;
    sector[i].site2ij = NULL;
    sector[i].ij2site = NULL;
    sector[i].propensity = NULL;
    sector[i].sites = NULL;
    sector[i].border = NULL;
  }

  // communicator needed between sweep sectors

  comm = new CommLattice2d(spk);
}

/* ---------------------------------------------------------------------- */

SweepLattice2d::~SweepLattice2d()
{
  delete random;
  delete comm;

  memory->destroy_2d_T_array(mask,nxlo,nylo);
  memory->destroy_2d_T_array(ranlat,1,1);

  for (int isector = 0; isector < nsector; isector++) {
    delete sector[isector].solve;
    memory->destroy_2d_T_array(sector[isector].site2ij);
    memory->destroy_2d_T_array(sector[isector].ij2site,nxlo,nylo);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->destroy_2d_T_array(sector[isector].border);
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice2d::init()
{
  int i,j,m;

  applattice = (AppLattice2d *) app;

  // error check

  if (Lmask && applattice->temperature > 0.0)
    error->all("Cannot mask sweeping with non-zero temperature");

  // get application params

  lattice = applattice->lattice;

  nx_local = applattice->nx_local;
  ny_local = applattice->ny_local;
  nx_offset = applattice->nx_offset;
  ny_offset = applattice->ny_offset;
  int nx_global = applattice->nx_global;
  int ny_global = applattice->ny_global;

  int procwest = applattice->procwest;
  int proceast = applattice->proceast;
  int procsouth = applattice->procsouth;
  int procnorth = applattice->procnorth;

  // app-specific settings

  delpropensity = applattice->delpropensity;
  delevent = applattice->delevent;
  numrandom = applattice->numrandom;

  nxlo = applattice->nxlo;
  nxhi = applattice->nxhi;
  nylo = applattice->nylo;
  nyhi = applattice->nyhi;

  delcol = delpropensity + 1;
  if (Lstrict) ncolor = delcol*delcol;
  else ncolor = 1;

  // setup sectors
  // partition sites by i,j values

  for (int isector = 0; isector < nsector; isector++) {
    delete sector[isector].solve;
    memory->destroy_2d_T_array(sector[isector].site2ij);
    memory->destroy_2d_T_array(sector[isector].ij2site);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->destroy_2d_T_array(sector[isector].border);
    sector[isector].solve = NULL;
    sector[isector].propensity = NULL;
    sector[isector].site2ij = NULL;
    sector[isector].ij2site = NULL;
    sector[isector].sites = NULL;
    sector[isector].border = NULL;
  }

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;

  sector[0].xlo = 1;
  sector[0].xhi = nx_half-1;
  sector[0].ylo = 1;
  sector[0].yhi = ny_half-1;
  sector[0].nx = sector[0].xhi - sector[0].xlo + 1;
  sector[0].ny = sector[0].yhi - sector[0].ylo + 1;

  sector[1].xlo = 1;
  sector[1].xhi = nx_half-1;
  sector[1].ylo = ny_half;
  sector[1].yhi = ny_local;
  sector[1].nx = sector[1].xhi - sector[1].xlo + 1;
  sector[1].ny = sector[1].yhi - sector[1].ylo + 1;
  
  sector[2].xlo = nx_half;
  sector[2].xhi = nx_local;
  sector[2].ylo = 1;
  sector[2].yhi = ny_half-1;
  sector[2].nx = sector[2].xhi - sector[2].xlo + 1;
  sector[2].ny = sector[2].yhi - sector[2].ylo + 1;

  sector[3].xlo = nx_half;
  sector[3].xhi = nx_local;
  sector[3].ylo = ny_half;
  sector[3].yhi = ny_local;
  sector[3].nx = sector[3].xhi - sector[3].xlo + 1;
  sector[3].ny = sector[3].yhi - sector[3].ylo + 1;

  // setup site2ij for each sector

  for (int isector = 0; isector < nsector; isector++) {
    int nsites = sector[isector].nx * sector[isector].ny;
    memory->create_2d_T_array(sector[isector].site2ij,nsites,2,
			      "sweep2d:site2ij");
    for (m = 0; m < nsites; m++) {
      i = m / sector[isector].ny + 1;
      j = m % sector[isector].ny +  1;
      sector[isector].site2ij[m][0] = i + sector[isector].xlo - 1;
      sector[isector].site2ij[m][1] = j + sector[isector].ylo - 1;
    }
  }

  // setup ij2site for each sector
  // 0 to nsite-1 for owned points in sector, else -1

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++) {
      memory->create_2d_T_array(sector[isector].ij2site,nxlo,nxhi,nylo,nyhi,
				"sweep2d:ij2site");
      for (i = 1-delpropensity; i <= nx_local+delpropensity; i++) 
	for (j = 1-delpropensity; j <= ny_local+delpropensity; j++)
	  sector[isector].ij2site[i][j] = -1;
      int nsites = sector[isector].nx * sector[isector].ny;
      for (m = 0; m < nsites; m++) {
	i = sector[isector].site2ij[m][0];
	j = sector[isector].site2ij[m][1];
	sector[isector].ij2site[i][j] = m;
      }
    }
  }

  // setup mask array for owned and ghost values
  // ghost values written to in app::site_event_rejection

  if (Lmask && mask == NULL) {
    memory->create_2d_T_array(mask,nxlo,nxhi,nylo,nyhi,"sweep2d:mask");
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= ny_local; j++)
	mask[i][j] = 0;
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
				"sweep2d:sites");
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_2d_T_array(ranlat,1,nx_local,1,ny_local,"sweep2d:ranlat");
    int isite;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++) {
	isite = (i+nx_offset)*ny_global + j + ny_offset;
	ranlat[i][j].init(seed+isite);
      }
  }

  // init communication for ghost sites

  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth,
	     delpropensity,delevent);

  // setup KMC solver and propensity arrays for each sector
  // propensity init requires all ghost cell info via comm->all()

  if (Lkmc) {
    comm->all(lattice);

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nx * sector[isector].ny;
      sector[isector].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep2d:propensity");
      for (m = 0; m < nsites; m++) {
	i = sector[isector].site2ij[m][0];
	j = sector[isector].site2ij[m][1];
	sector[isector].propensity[m] = applattice->site_propensity(i,j);
      }

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }

  // If adaptive kmc specified
  // compute deln0, which controls future timesteps
  // If deln0 arleady specified, compute delt

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
    stoptime = applattice->stoptime;
    delt = MIN(delt,stoptime);
  }
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice2d::do_sweep(double &dt)
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

  if (Lkmc && Ladapt) {
    MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_SUM,world);
    if (pmaxall > 0.0) delt = deln0/pmaxall;
    delt = MIN(delt,stoptime);
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector(int icolor, int isector)
{
  int i,j;
  
  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  
  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      applattice->site_event_rejection(i,j,random);
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites with masking
   perform events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_mask(int icolor, int isector)
{
  int i,j;
  
  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  
  boundary_clear_mask(isector);
  
  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++) {
      if (mask[i][j]) continue;
      applattice->site_event_rejection(i,j,random);
    }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_strict(int icolor, int isector)
{
  int i,j,i0,j0;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  i0 = icolor/delcol  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor%delcol  - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  for (i = xlo+i0; i <= xhi; i += delcol)
    for (j = ylo+j0; j <= yhi; j += delcol)
      applattice->site_event_rejection(i,j,&ranlat[i][j]);
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color with masking
   perform events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_mask_strict(int icolor, int isector)
{
  int i,j,m,i0,j0;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  boundary_clear_mask(isector);
  
  i0 = icolor/delcol  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor%delcol  - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  for (i = xlo+i0; i <= xhi; i += delcol)
    for (j = ylo+j0; j <= yhi; j += delcol) {
      if (mask[i][j]) {
	for (m = 0; m < numrandom; m++) ranlat[i][j].uniform();
	continue;
      }
      applattice->site_event_rejection(i,j,&ranlat[i][j]);
    }
}

/* ----------------------------------------------------------------------
   generate events in each sector using KMC solver
 ------------------------------------------------------------------------- */

void SweepLattice2d::sweep_sector_kmc(int icolor, int isector)
{
  double dt,time;
  int done,isite,i,j;

  // extract sector specific info from quad struct

  int **border = sector[isector].border;
  int nborder = sector[isector].nborder;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int **site2ij = sector[isector].site2ij;
  int **ij2site = sector[isector].ij2site;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // solver, propensity, ij2site array

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  applattice->propensity = propensity;
  applattice->ij2site = ij2site;

  // update propensities for sites which neighbor a site outside sector
  // necessary since outside sites may have changed
  // attribute this chunk of time to comm, b/c due to decomposition

  timer->stamp();

  int nsites = 0;

  for (int m = 0; m < nborder; m++) {
    i = border[m][0];
    j = border[m][1];
    isite = ij2site[i][j];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j);
  }
    
  solve->update(nsites,sites,propensity);
  timer->stamp(TIME_COMM);

  // compute maximum sector propensity per site
  // controls future timesteps. We do this at start of sweep
  // for the simple reason that a sector propensity can drop to
  // zero during a sweep, but it can never increase from zero.
  // Hence we avoid the problem of all the sectors reporting
  // zero propensity. a more correct way to handle this would be to do
  // a full propensity update when all the sectors have been swept,
  // but that would increase the cost substantially.
  
  if (Ladapt) {
    int ntmp = applattice->solve->get_num_active();
    if (ntmp > 0) {
      double ptmp = applattice->solve->get_total_propensity();
      ptmp /= ntmp;
      pmax = MAX(ptmp,pmax);
    }
  }

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
	i = site2ij[isite][0];
	j = site2ij[isite][1];
	applattice->site_event(i,j,0,random);
	applattice->ntimestep++;
      }
      timer->stamp(TIME_APP);
    }
  }
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = NULL;
  applattice->ij2site = NULL;
}

/* ----------------------------------------------------------------------
   create list of border sites for a sector
   border site = site in sector with a neighbor outside the sector
   neighbor can be another owned site (outside sector) or a ghost
   border = lattice indices of the sites
 ------------------------------------------------------------------------- */

int SweepLattice2d::find_border_sites(int isector)
{
  int **border = NULL;
  int nborder = 0;
  int nmax = 0;

  // loop over all sites in sector
  // only add those within delta of boundary to border

  int i,j;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  int delta = delpropensity + delevent;

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++) {
      if (i-xlo >= delta && xhi-i >= delta && 
	  j-ylo >= delta && yhi-j >= delta) continue;

      if (nborder == nmax) {
	nmax += DELTA;
	memory->grow_2d_T_array(border,nmax,2,"sweep2d:border");
      }
      border[nborder][0] = i;
      border[nborder][1] = j;
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

void SweepLattice2d::boundary_clear_mask(int isector)
{
  int **border = sector[isector].border;
  int nborder = sector[isector].nborder;

  for (int m = 0; m < nborder; m++)
    mask[border[m][0]][border[m][1]] = 0;
}
