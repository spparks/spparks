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
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SweepLattice3d::SweepLattice3d(SPK *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
  
  // parse optional args

  Lmask = true;
  mask = NULL;
  Lstrict = false;
  Lpicklocal = false;
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
    } else if (strcmp(arg[iarg],"delt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      delt = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal sweep_style command");
  }

  // set sweep method function ptr

  if (Lstrict) {
    if (Lmask) {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sector = &SweepLattice3d::sweep_quadrant_mask_strict;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sector = &SweepLattice3d::sweep_quadrant_strict;
    }
  } else {
    if (Lmask) {
      if (Lpicklocal) sector = &SweepLattice3d::sweep_quadrant_mask_picklocal;
      else sector = &SweepLattice3d::sweep_quadrant_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sector = &SweepLattice3d::sweep_quadrant;
    }
  }

  if (Lstrict) ncolor = 8;
  else ncolor = 1;

  // communicator needed between sweep sectors

  comm = new CommLattice3d(spk);
}

/* ---------------------------------------------------------------------- */

SweepLattice3d::~SweepLattice3d()
{
  delete random;
  delete comm;

  memory->destroy_3d_T_array(mask);
  memory->destroy_3d_T_array(ranlat);
}

/* ---------------------------------------------------------------------- */

void SweepLattice3d::init()
{
  applattice = (AppLattice3d *) app;

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

  temperature = applattice->temperature;
  if (temperature != 0.0) t_inverse = 1.0/temperature;
  masklimit = applattice->masklimit;
  
  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

  nquad = 8;

  quad[0].xlo = 1;
  quad[0].xhi = nx_half-1;
  quad[0].ylo = 1;
  quad[0].yhi = ny_half-1;
  quad[0].zlo = 1;
  quad[0].zhi = nz_half-1;
  
  quad[1].xlo = 1;
  quad[1].xhi = nx_half-1;
  quad[1].ylo = 1;
  quad[1].yhi = ny_half-1;
  quad[1].zlo = nz_half;
  quad[1].zhi = nz_local;

  quad[2].xlo = 1;
  quad[2].xhi = nx_half-1;
  quad[2].ylo = ny_half;
  quad[2].yhi = ny_local;
  quad[2].zlo = 1;
  quad[2].zhi = nz_half-1;

  quad[3].xlo = 1;
  quad[3].xhi = nx_half-1;
  quad[3].ylo = ny_half;
  quad[3].yhi = ny_local;
  quad[3].zlo = nz_half;
  quad[3].zhi = nz_local;

  quad[4].xlo = nx_half;
  quad[4].xhi = nx_local;
  quad[4].ylo = 1;
  quad[4].yhi = ny_half-1;
  quad[4].zlo = 1;
  quad[4].zhi = nz_half-1;

  quad[5].xlo = nx_half;
  quad[5].xhi = nx_local;
  quad[5].ylo = 1;
  quad[5].yhi = ny_half-1;
  quad[5].zlo = nz_half;
  quad[5].zhi = nz_local;

  quad[6].xlo = nx_half;
  quad[6].xhi = nx_local;
  quad[6].ylo = ny_half;
  quad[6].yhi = ny_local;
  quad[6].zlo = 1;
  quad[6].zhi = nz_half-1;

  quad[7].xlo = nx_half;
  quad[7].xhi = nx_local;
  quad[7].ylo = ny_half;
  quad[7].yhi = ny_local;
  quad[7].zlo = nz_half;
  quad[7].zhi = nz_local;

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    memory->create_3d_T_array(mask,nx_local+2,ny_local+2,nz_local+2,
			      "sweeplattice3d:mask");
    
    for (int i = 1; i <= nx_local; i++) 
      for (int j = 1; j <= ny_local; j++) 
	for (int k = 1; k <= nz_local; k++) 
	  mask[i][j][k] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_3d_T_array(ranlat,nx_local+2,ny_local+2,nz_local+2,
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

  // init communication for ghost sites
  
  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup);
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice3d::do_sweep(double &dt)
{
  for (int icolor = 0; icolor < ncolor; icolor++)
    for (int iquad = 0; iquad < nquad; iquad++) {
      timer->stamp();
      comm->sector(lattice,iquad);
      timer->stamp(TIME_COMM);

      (this->*sector)(icolor,iquad);
      timer->stamp(TIME_SOLVE);
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
   
void SweepLattice3d::sweep_quadrant(int icolor, int iquad)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

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
   
void SweepLattice3d::sweep_quadrant_mask(int icolor, int iquad)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

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
	
	if (efinal <= einitial) continue;
	else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
      }
}

/* ---------------------------------------------------------------------- */

void SweepLattice3d::sweep_quadrant_mask_picklocal(int icolor, int iquad)
{
  int i,j,k,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

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
   
void SweepLattice3d::sweep_quadrant_strict(int icolor, int iquad)
{
  int i,j,k,i0,j0,k0,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  i0 = (icolor/4 + nx_offset + xlo) % 2;
  j0 = (icolor/2 + ny_offset + ylo) % 2;
  k0 = (icolor   + nz_offset + zlo) % 2;

  for (i = xlo+i0; i <= xhi; i += 2)
    for (j = ylo+j0; j <= yhi; j += 2)
      for (k = zlo+k0;  k <= zhi; k+=2) {
	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);
	
	newstate = 
	  applattice->site_pick_random(i,j,k,ranlat[i][j][k].uniform());
	lattice[i][j][k] = newstate;
	efinal = applattice->site_energy(i,j,k);
	
	if (efinal <= einitial) continue;
	else if (temperature == 0.0) lattice[i][j][k] = oldstate;
	else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	  lattice[i][j][k] = oldstate;
      }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice3d::sweep_quadrant_mask_strict(int icolor, int iquad)
{
  int i,j,k,i0,j0,k0,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

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

  i0 = (icolor/4 + nx_offset + xlo) % 2;
  j0 = (icolor/2 + ny_offset + ylo) % 2;
  k0 = (icolor   + nz_offset + zlo) % 2;

  // call RNG even if skipping via mask
  // insures same answer no matter where proc boundaries are

  for (i = xlo+i0; i <= xhi; i+=2)
    for (j = ylo+j0; j <= yhi; j+=2)
      for (k = zlo+k0;  k <= zhi; k+=2) {
	if (mask[i][j][k]) {
	  ranlat[i][j][k].uniform();
	  continue;
	}

	oldstate = lattice[i][j][k];
	einitial = applattice->site_energy(i,j,k);
	if (einitial < masklimit) {
	  mask[i][j][k] = 1;
	  ranlat[i][j][k].uniform();
	  continue;
	}
	
	newstate =
	  applattice->site_pick_random(i,j,k,ranlat[i][j][k].uniform());
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
