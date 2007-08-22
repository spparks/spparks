/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SweepLattice2d::SweepLattice2d(SPK *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  
  // parse optional args

  Lmask = false;
  mask = NULL;
  Lstrict = false;
  Lpicklocal = false;
  Lkmc = false;
  ranlat = NULL;
  delt = 1.0;

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
    } else if (strcmp(arg[iarg],"picklocal") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lpicklocal = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lpicklocal = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kmc") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lkmc = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lkmc = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      delt = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal sweep_style command");
  }

  // set sweep method function ptr

  if (Lkmc) {
    if (Lmask || Lpicklocal)
      error->all("Combination of sweep flags is unsupported");
    sweeper = &SweepLattice2d::sweep_sector_kmc;
  } else if (Lstrict) {
    if (Lmask) {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice2d::sweep_sector_mask_strict;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice2d::sweep_sector_strict;
    }
  } else {
    if (Lmask) {
      if (Lpicklocal) sweeper = &SweepLattice2d::sweep_sector_mask_picklocal;
      else sweeper = &SweepLattice2d::sweep_sector_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice2d::sweep_sector;
    }
  }

  // initializations

  nsector = 4;
  for (int i = 0; i < nsector; i++) {
    sector[i].solve = NULL;
    sector[i].propensity = NULL;
    sector[i].site2ij = NULL;
    sector[i].sites = NULL;
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

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++) {
      delete sector[isector].solve;
      memory->sfree(sector[isector].propensity);
      memory->destroy_2d_T_array(sector[isector].site2ij);
      memory->sfree(sector[isector].sites);
    }
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice2d::init()
{
  applattice = (AppLattice2d *) app;

  lattice = applattice->lattice;
  ij2site = applattice->ij2site;

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

  temperature = applattice->temperature;
  if (temperature != 0.0) t_inverse = 1.0/temperature;

  // app-specific settings

  masklimit = applattice->masklimit;
  delghost = applattice->delghost;
  dellocal = applattice->dellocal;
  nxlo = applattice->nxlo;
  nxhi = applattice->nxhi;
  nylo = applattice->nylo;
  nyhi = applattice->nyhi;

  delcol = delghost+1;
  if (Lstrict) ncolor = delcol*delcol;
  else ncolor = 1;

  // setup sectors

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;

  sector[0].xlo = 1;
  sector[0].xhi = nx_half-1;
  sector[0].ylo = 1;
  sector[0].yhi = ny_half-1;
  sector[0].nx = sector[0].xhi - sector[0].xlo + 1;
  sector[0].ny = sector[0].yhi - sector[0].ylo + 1;
  delete sector[0].solve;
  memory->sfree(sector[0].propensity);
  memory->destroy_2d_T_array(sector[0].site2ij);
  memory->sfree(sector[0].sites);

  sector[1].xlo = 1;
  sector[1].xhi = nx_half-1;
  sector[1].ylo = ny_half;
  sector[1].yhi = ny_local;
  sector[1].nx = sector[1].xhi - sector[1].xlo + 1;
  sector[1].ny = sector[1].yhi - sector[1].ylo + 1;
  delete sector[1].solve;
  memory->sfree(sector[1].propensity);
  memory->destroy_2d_T_array(sector[1].site2ij);
  memory->sfree(sector[1].sites);
  
  sector[2].xlo = nx_half;
  sector[2].xhi = nx_local;
  sector[2].ylo = 1;
  sector[2].yhi = ny_half-1;
  sector[2].nx = sector[2].xhi - sector[2].xlo + 1;
  sector[2].ny = sector[2].yhi - sector[2].ylo + 1;
  delete sector[2].solve;
  memory->sfree(sector[2].propensity);
  memory->destroy_2d_T_array(sector[2].site2ij);
  memory->sfree(sector[2].sites);

  sector[3].xlo = nx_half;
  sector[3].xhi = nx_local;
  sector[3].ylo = ny_half;
  sector[3].yhi = ny_local;
  sector[3].nx = sector[3].xhi - sector[3].xlo + 1;
  sector[3].ny = sector[3].yhi - sector[3].ylo + 1;
  delete sector[3].solve;
  memory->sfree(sector[3].propensity);
  memory->destroy_2d_T_array(sector[3].site2ij);
  memory->sfree(sector[3].sites);

  // init communication for ghost sites

  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth,
	     delghost,dellocal);

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    memory->create_2d_T_array(mask,nxlo,nxhi,nylo,nyhi,
			      "sweeplattice2d:mask");
    for (int i = 1; i <= nx_local; i++) 
      for (int j = 1; j <= ny_local; j++) 
	mask[i][j] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_2d_T_array(ranlat,1,nx_local,1,ny_local,
			      "sweeplattice2d:ranlat");
    int isite;
    for (int i = 1; i <= nx_local; i++)
      for (int j = 1; j <= ny_local; j++) {
	isite = (i+nx_offset)*ny_global + j + ny_offset;
	ranlat[i][j].init(seed+isite);
      }
  }

  // setup KMC solver and propensity arrays for each sector
  // set ij2site and site2ij values in app to reflect sector mapping
  // propensity init requires ghost cell info for entire sub-domain

  if (Lkmc) {
    int i,j,m;

    comm->all(lattice);

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nx * sector[isector].ny;
      int nborder = 2*sector[isector].nx + 2*sector[isector].ny;
      sector[isector].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep:propensity");
      memory->create_2d_T_array(sector[isector].site2ij,nsites,2,
				"sweep:site2ij");
      sector[isector].sites =
	(int*) memory->smalloc(nborder*sizeof(int),"sweep:sites");

      for (i = sector[isector].xlo; i <= sector[isector].xhi; i++)
	for (j = sector[isector].ylo; j <= sector[isector].yhi; j++)
	  ij2site[i][j] = 
	    (i-sector[isector].xlo)*sector[isector].ny + j-sector[isector].ylo;

      for (m = 0; m < nsites; m++) {
	i = m / sector[isector].ny + 1;
	j = m % sector[isector].ny +  1;
	sector[isector].site2ij[m][0] = i + sector[isector].xlo - 1;
	sector[isector].site2ij[m][1] = j + sector[isector].ylo - 1;
      }

      for (i = sector[isector].xlo; i <= sector[isector].xhi; i++)
	for (j = sector[isector].ylo; j <= sector[isector].yhi; j++)
	  sector[isector].propensity[ij2site[i][j]] = 
	    applattice->site_propensity(i,j,0);

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice2d::do_sweep(double &dt)
{
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
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   application picks a new state for the site
   compute energy change due to state change
   no energy change = accept
   downhill energy change = accept
   uphill energy change = accept/reject via Boltzmann factor
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector(int icolor, int isector)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++) {
      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);

      newstate = applattice->site_pick_random(i,j,random->uniform());
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) continue;
      else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
    }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   skip sites that can't change via mask
 ------------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_mask(int icolor, int isector)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi);
  
  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++) {
      if (mask[i][j]) continue;

      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);
      if (einitial < masklimit) {
	mask[i][j] = 1;
	continue;
      }

      newstate = applattice->site_pick_random(i,j,random->uniform());
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) {
	applattice->site_clear_mask(mask,i,j);
	continue;
      } else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
      else applattice->site_clear_mask(mask,i,j);
    }
}

/* ---------------------------------------------------------------------- */

void SweepLattice2d::sweep_sector_mask_picklocal(int icolor, int isector)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi);
  
  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++) {
      if (mask[i][j]) continue;

      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);
      if (einitial < masklimit) {
	mask[i][j] = 1;
	continue;
      }

      newstate = applattice->site_pick_local(i,j,random->uniform());
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) {
	applattice->site_clear_mask(mask,i,j);
	continue;
      } else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
      else applattice->site_clear_mask(mask,i,j);
    }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_strict(int icolor, int isector)
{
  int i,j,i0,j0,oldstate,newstate;
  double einitial,efinal;
  double xtmp1,xtmp2;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  i0 = icolor/delcol  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor%delcol  - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  for (i = xlo+i0; i <= xhi; i += delcol)
    for (j = ylo+j0; j <= yhi; j += delcol) {

      // call RNG even if skipping via mask
      // insures same answer no matter where proc boundaries are

      xtmp1 = ranlat[i][j].uniform();
      xtmp2 = ranlat[i][j].uniform();

      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);

      newstate = applattice->site_pick_random(i,j,xtmp1);
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) continue;

      else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (xtmp2 > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
    }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_sector_mask_strict(int icolor, int isector)
{
  int i,j,i0,j0,oldstate,newstate;
  double einitial,efinal;
  double xtmp1,xtmp2;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi);
  
  i0 = icolor/delcol  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor%delcol  - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  for (i = xlo+i0; i <= xhi; i+=delcol)
    for (j = ylo+j0; j <= yhi; j+=delcol) {

  // call RNG even if skipping via mask
  // insures same answer no matter where proc boundaries are

      xtmp1 = ranlat[i][j].uniform();
      xtmp2 = ranlat[i][j].uniform();

      if (mask[i][j]) continue;

      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);
      if (einitial < masklimit) {
	mask[i][j] = 1;
	continue;
      }

      newstate = applattice->site_pick_random(i,j,xtmp1);
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) {
	applattice->site_clear_mask(mask,i,j);
	continue;
      } else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (xtmp2 > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
      else applattice->site_clear_mask(mask,i,j);
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

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int **site2ij = sector[isector].site2ij;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // sector bounds, propensity array, solver

  applattice->nx_sector_lo = xlo;
  applattice->nx_sector_hi = xhi;
  applattice->ny_sector_lo = ylo;
  applattice->ny_sector_hi = yhi;

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  double *hold_propensity = applattice->propensity;
  applattice->propensity = propensity;

  // update propensities on all 4 sector edges
  // necessary since ghosts of sector may have changed

  int nsites = 0;

  j = ylo;
  for (i = xlo; i <= xhi; i++) {
    isite = ij2site[i][j];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j,0);
  }
  j = yhi;
  for (i = xlo; i <= xhi; i++) {
    isite = ij2site[i][j];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j,0);
  }

  i = xlo;
  for (j = ylo; j <= yhi; j++) {
    isite = ij2site[i][j];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j,0);
  }
  i = xhi;
  for (j = ylo; j <= yhi; j++) {
    isite = ij2site[i][j];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,j,0);
  }

  solve->update(nsites,sites,propensity);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    // do not allow threshold time to be exceeded
    time += dt;	
    if (isite < 0 || time >= delt) done = 1;
    else {
      i = site2ij[isite][0];
      j = site2ij[isite][1];
      applattice->site_event(i,j,0);
      timer->stamp(TIME_APP);
    }
  }
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = hold_propensity;
}

/* ----------------------------------------------------------------------
   unset all masks on boundary
   may be out of date, due to state change on neighboring processor
   could reverse comm mask values, but that might be slower
 ------------------------------------------------------------------------- */

void SweepLattice2d::boundary_clear_mask(int xlo, int xhi, int ylo, int yhi) {
  int i,j;

  if (delghost == 1) {
    if (ylo == 1) {
      for (i = xlo; i <= xhi; i++)
	mask[i][1] = 0;
    } else if (yhi == ny_local)
      for (i = xlo; i <= xhi; i++)
	mask[i][ny_local] = 0;
    
    if (xlo == 1) {
      for (j = ylo; j <= yhi; j++)
	mask[1][j] = 0;
    } else if (xhi == nx_local)
      for (j = ylo; j <= yhi; j++)
	mask[nx_local][j] = 0;
  } else {
    if (ylo == 1) {
      for (i = xlo; i <= xhi; i++)
	for (j = ylo; j < ylo+delghost; j++)
	  mask[i][j] = 0;
    } else if (yhi == ny_local)
      for (i = xlo; i <= xhi; i++)
	for (j = yhi; j > yhi-delghost; j--)
	  mask[i][j] = 0;
    
    if (xlo == 1) {
      for (i = xlo; i < xlo+delghost; i++)
	for (j = ylo; j <= yhi; j++)
	  mask[i][j] = 0;
    } else if (xhi == nx_local)
      for (i = xhi; i > xhi-delghost; i--)
	for (j = ylo; j <= yhi; j++)
	  mask[i][j] = 0;
  }
}
