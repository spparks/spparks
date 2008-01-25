/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SweepLattice3d::SweepLattice3d(SPK *spk, int narg, char **arg) : 
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
    sweeper = &SweepLattice3d::sweep_sector_kmc;
  } else if (Lstrict) {
    if (Lmask) {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice3d::sweep_sector_mask_strict;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice3d::sweep_sector_strict;
    }
  } else {
    if (Lmask) {
      if (Lpicklocal) sweeper = &SweepLattice3d::sweep_sector_mask_picklocal;
      else sweeper = &SweepLattice3d::sweep_sector_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice3d::sweep_sector;
    }
  }

  nsector = 8;
  for (int i = 0; i < nsector; i++) {
    sector[i].solve = NULL;
    sector[i].propensity = NULL;
    sector[i].site2ijk = NULL;
    sector[i].sites = NULL;
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
    memory->sfree(sector[isector].propensity);
    memory->destroy_2d_T_array(sector[isector].site2ijk);
    memory->sfree(sector[isector].sites);
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice3d::init()
{
  applattice = (AppLattice3d *) app;

  lattice = applattice->lattice;
  ijk2site = applattice->ijk2site;

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
  nzlo = applattice->nzlo;
  nzhi = applattice->nzhi;

  delcol = delghost+1;
  if (Lstrict) ncolor = delcol*delcol*delcol;
  else ncolor = 1;
  
  // setup sectors

  for (int i = 0; i < nsector; i++) {
    delete sector[i].solve;
    memory->sfree(sector[i].propensity);
    memory->destroy_2d_T_array(sector[i].site2ijk);
    memory->sfree(sector[i].sites);
    sector[i].solve = NULL;
    sector[i].propensity = NULL;
    sector[i].site2ijk = NULL;
    sector[i].sites = NULL;
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

  // init communication for ghost sites

  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup,
	     delghost,dellocal);

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    memory->create_3d_T_array(mask,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
			      "sweeplattice3d:mask");
    
    for (int i = 1; i <= nx_local; i++) 
      for (int j = 1; j <= ny_local; j++) 
	for (int k = 1; k <= nz_local; k++) 
	  mask[i][j][k] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_3d_T_array(ranlat,1,nx_local,1,ny_local,1,nz_local,
			      "sweeplattice3d:ranlat");
    int isite;
    for (int i = 1; i <= nx_local; i++)
      for (int j = 1; j <= ny_local; j++)
	for (int k = 1; k <= nz_local; k++) {
	  isite = (i+nx_offset)*ny_global*nz_global + 
	    (j+ny_offset)*nz_global + k + nz_offset;
	  ranlat[i][j][k].init(seed+isite);
	}
  }

  // setup KMC solver and propensity arrays for each sector
  // reset app's ij2site values to reflect sector mapping
  // create site2ij values for each sector
  // propensity init requires ghost cell info for entire sub-domain

  if (Lkmc) {
    int i,j,k,m;

    psum = 0.0;

    comm->all(lattice);

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nx * sector[isector].ny * 
	sector[isector].nz;
      sector[isector].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep:propensity");
      memory->create_2d_T_array(sector[isector].site2ijk,nsites,3,
				"sweep:ijk2site");

      for (i = sector[isector].xlo; i <= sector[isector].xhi; i++)
	for (j = sector[isector].ylo; j <= sector[isector].yhi; j++)
	  for (k = sector[isector].zlo; k <= sector[isector].zhi; k++)
	    ijk2site[i][j][k] = 
	      (i-sector[isector].xlo)*sector[isector].ny*sector[isector].nz +
	      (j-sector[isector].ylo)*sector[isector].nz + 
	      k-sector[isector].zlo;

      int nborder = 2*(dellocal+delghost) * 
	(sector[isector].nx*sector[isector].ny + 
	 sector[isector].ny*sector[isector].nz + 
	 sector[isector].nx*sector[isector].nz);
      sector[isector].sites =
	(int*) memory->smalloc(nborder*sizeof(int),"sweep:sites");

      for (m = 0; m < nsites; m++) {
	i = m / sector[isector].ny/sector[isector].nz + 1;
	j = (m / sector[isector].nz) % sector[isector].ny + 1;
	k = m % sector[isector].nz + 1;
	sector[isector].site2ijk[m][0] = i + sector[isector].xlo - 1;
	sector[isector].site2ijk[m][1] = j + sector[isector].ylo - 1;
	sector[isector].site2ijk[m][2] = k + sector[isector].zlo - 1;
      }

      for (i = sector[isector].xlo; i <= sector[isector].xhi; i++)
	for (j = sector[isector].ylo; j <= sector[isector].yhi; j++)
	  for (k = sector[isector].zlo; k <= sector[isector].zhi; k++)
	    sector[isector].propensity[ijk2site[i][j][k]] = 
	      applattice->site_propensity(i,j,k,0);

      sector[isector].solve->init(nsites,sector[isector].propensity);
      psum += sector[isector].solve->get_total_propensity();
    }

    // Compute deln0, which controls future timesteps
    MPI_Allreduce(&psum,&psumall,1,MPI_DOUBLE,MPI_SUM,world);
    deln0 = psumall*delt;
  }

}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice3d::do_sweep(double &dt)
{
  if (Lkmc) {
    psum = 0.0;
  }

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
  if (Lkmc) {
    MPI_Allreduce(&psum,&psumall,1,MPI_DOUBLE,MPI_SUM,world);
    if (psumall > 0.0) {
      delt = deln0/psumall;
    }
  }

}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   application picks a new state for the site
   compute energy change due to state change
   no energy change = accept
   downhill energy change = accept
   uphill energy change = accept/reject via Boltzmann factor
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector(int icolor, int isector)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);

	newstate = applattice->site_pick_random(i,j,k,random->uniform());
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);

	if (efinal <= einitial) continue;
	else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
      }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   skip sites that can't change via mask
 ------------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_mask(int icolor, int isector)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi,zlo,zhi);

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	if (mask[i][j][k]) continue;

	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);

	if (einitial < masklimit) {
	  mask[i][j][k] = 1;
	  continue;
	}
	
	newstate = applattice->site_pick_random(i,j,k,random->uniform());
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);
	
	if (efinal <= einitial) {
	  applattice->site_clear_mask(mask,i,j,k);
	  continue;
	} else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
	else applattice->site_clear_mask(mask,i,j,k);
      }
}

/* ---------------------------------------------------------------------- */

void SweepLattice3d::sweep_sector_mask_picklocal(int icolor, int isector)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi,zlo,zhi);

  for (i = xlo; i <= xhi; i++)
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	if (mask[i][j][k]) continue;

	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);
	if (einitial < masklimit) {
	  mask[i][j][k] = 1;
	  continue;
	}
	
	newstate = applattice->site_pick_local(i,j,k,random->uniform());
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);
	
	if (efinal <= einitial) {
	  applattice->site_clear_mask(mask,i,j,k);
	  continue;
	} else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
	else applattice->site_clear_mask(mask,i,j,k);
      }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_strict(int icolor, int isector)
{
  int i,j,k,i0,j0,k0,oldstate,newstate;
  double einitial,efinal;
  double xtmp1,xtmp2;

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
      for (k = zlo+k0;  k <= zhi; k += delcol) {

	// call RNG even if skipping via mask
	// insures same answer no matter where proc boundaries are
	
	xtmp1 = ranlat[i][j][k].uniform();
	xtmp2 = ranlat[i][j][k].uniform();

	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);
	
	newstate = 
	  applattice->site_pick_random(i,j,k,xtmp1);
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);
	
	if (efinal <= einitial) continue;
	else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (xtmp2 > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
      }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_sector_mask_strict(int icolor, int isector)
{
  int i,j,k,i0,j0,k0,oldstate,newstate;
  double einitial,efinal;
  double xtmp1,xtmp2;

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  boundary_clear_mask(xlo,xhi,ylo,yhi,zlo,zhi);

  i0 = icolor/(delcol*delcol)  - (nx_offset + xlo-1) % delcol;
  i0 = (i0 < 0) ? i0+delcol : i0;

  j0 = icolor/delcol % delcol - (ny_offset + ylo-1) % delcol;
  j0 = (j0 < 0) ? j0+delcol : j0;

  k0 = icolor%delcol  - (nz_offset + zlo-1) % delcol;
  k0 = (k0 < 0) ? k0+delcol : k0;

  for (i = xlo+i0; i <= xhi; i+=delcol)
    for (j = ylo+j0; j <= yhi; j+=delcol)
      for (k = zlo+k0;  k <= zhi; k += delcol) {

	// call RNG even if skipping via mask
	// insures same answer no matter where proc boundaries are
	
	xtmp1 = ranlat[i][j][k].uniform();
	xtmp2 = ranlat[i][j][k].uniform();

	if (mask[i][j][k]) continue;

	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);
	if (einitial < masklimit) {
	  mask[i][j][k] = 1;
	  continue;
	}
	
	newstate =
	  applattice->site_pick_random(i,j,k,xtmp1);
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);
	
	if (efinal <= einitial) {
	  applattice->site_clear_mask(mask,i,j,k);
	  continue;
	} else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (xtmp2 > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
	else applattice->site_clear_mask(mask,i,j,k);
      }
}

/* ----------------------------------------------------------------------
   generate events in each sector using KMC solver
 ------------------------------------------------------------------------- */

void SweepLattice3d::sweep_sector_kmc(int icolor, int isector)
{
  double dt,time;
  int done,isite,i,j,k;

  timer->stamp();

  // extract sector specific info from octant struct

  int xlo = sector[isector].xlo;
  int xhi = sector[isector].xhi;
  int ylo = sector[isector].ylo;
  int yhi = sector[isector].yhi;
  int zlo = sector[isector].zlo;
  int zhi = sector[isector].zhi;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int **site2ijk = sector[isector].site2ijk;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // sector bounds, solver, propensity array

  applattice->nx_sector_lo = xlo;
  applattice->nx_sector_hi = xhi;
  applattice->ny_sector_lo = ylo;
  applattice->ny_sector_hi = yhi;
  applattice->nz_sector_lo = zlo;
  applattice->nz_sector_hi = zhi;

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  double *hold_propensity = applattice->propensity;
  applattice->propensity = propensity;

  // update owned propensities on all 6 sector faces
  // necessary since sector ghosts may have changed

  int nsites = 0;
  int deltemp = dellocal+delghost;

  if (deltemp <= 1) {
    k = zlo;
    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= yhi; j++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }
    k = zhi;
    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= yhi; j++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }

    j = ylo;
    for (i = xlo; i <= xhi; i++)
      for (k = zlo; k <= zhi; k++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }
    j = yhi;
    for (i = xlo; i <= xhi; i++)
      for (k = zlo; k <= zhi; k++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }
    
    i = xlo;
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }
    i = xhi;
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++) {
	isite = ijk2site[i][j][k];
	sites[nsites++] = isite;
	propensity[isite] = applattice->site_propensity(i,j,k,0);
      }

  } else {
    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= yhi; j++)
	for (k = zlo; k <= zlo+deltemp-1; k++) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}
    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= yhi; j++)
	for (k = zhi; k >= zhi-deltemp+1; k--) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}

    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= ylo+deltemp-1; j++)
	for (k = zlo; k <= zhi; k++) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}
    for (i = xlo; i <= xhi; i++)
      for (j = yhi; j >= yhi-deltemp+1; j--)
	for (k = zlo; k <= zhi; k++) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}
    
    for (i = xlo; i <= xlo+deltemp-1; i++)
      for (j = ylo; j <= yhi; j++)
	for (k = zlo; k <= zhi; k++) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}
    for (i = xhi; i >= xhi-deltemp+1; i--)
      for (j = ylo; j <= yhi; j++)
	for (k = zlo; k <= zhi; k++) {
	  isite = ijk2site[i][j][k];
	  sites[nsites++] = isite;
	  propensity[isite] = applattice->site_propensity(i,j,k,0);
	}
  }

  solve->update(nsites,sites,propensity);

  // attribute this chunk of time to Communication,
  // because it is due to the spatial decomposition
  timer->stamp(TIME_COMM);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
    applattice->ntimestep++;
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);
    if (isite < 0) done = 1;
    else {
      i = site2ijk[isite][0];
      j = site2ijk[isite][1];
      k = site2ijk[isite][2];
      applattice->site_event(i,j,k,0);
      time += dt;	
      if (time >= delt) done = 1;
      timer->stamp(TIME_APP);
    }
  }

  // Sum the final propensity

  psum += applattice->solve->get_total_propensity();
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = hold_propensity;
}

/* ----------------------------------------------------------------------
   unset all masks on boundary
   may be out of date, due to state change on neighboring processor
   could reverse comm mask values, but that might be slower
 ------------------------------------------------------------------------- */

void SweepLattice3d::boundary_clear_mask(int xlo, int xhi, int ylo, 
					 int yhi, int zlo, int zhi) {
  int i,j,k;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

  if (delghost == 1) {
    if (zlo == 1) k = zlo;
    else if (zhi == nz_local) k = zhi;
    for (i = xlo; i <= xhi; i++)
      for (j = ylo; j <= yhi; j++)
	mask[i][j][k] = 0;

    if (ylo == 1) j = ylo;
    else if (yhi == ny_local) j = yhi;
    for (i = xlo; i <= xhi; i++)
      for (k = zlo; k <= zhi; k++)
	mask[i][j][k] = 0;

    if (xlo == 1) i = xlo;
    else if (xhi == nx_local) i = xhi;
    for (j = ylo; j <= yhi; j++)
      for (k = zlo; k <= zhi; k++)
	mask[i][j][k] = 0;
  } else {
    if (zlo == 1)
      for (i = xlo; i <= xhi; i++)
	for (j = ylo; j <= yhi; j++)
	  for (k = zlo; k < zlo+delghost; k++)
	    mask[i][j][k] = 0;
    else if (zhi == nz_local) 
      for (i = xlo; i <= xhi; i++)
	for (j = ylo; j <= yhi; j++)
	  for (k = zhi; k > zlo-delghost; k--)
	    mask[i][j][k] = 0;

    if (ylo == 1)
      for (i = xlo; i <= xhi; i++)
	for (j = ylo; j < ylo+delghost; j++)
	  for (k = zlo; k <= zhi; k++)
	    mask[i][j][k] = 0;
    else if (yhi == nz_local) 
      for (i = xlo; i <= xhi; i++)
	for (j = yhi; j > ylo-delghost; j--)
	  for (k = zlo; k <= zhi; k++)
	    mask[i][j][k] = 0;

    if (xlo == 1)
      for (i = xlo; i < xlo+delghost; i++)
	for (j = ylo; j <= yhi; j++)
	  for (k = zlo; k <= zhi; k++)
	    mask[i][j][k] = 0;
    else if (xhi == nx_local) 
      for (i = xhi; i > xlo-delghost; i--)
	for (j = ylo; j <= yhi; j++)
	  for (k = zlo; k <= zhi; k++)
	    mask[i][j][k] = 0;

  }
}
