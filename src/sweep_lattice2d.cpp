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
      else sector = &SweepLattice2d::sweep_quadrant_mask_strict;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sector = &SweepLattice2d::sweep_quadrant_strict;
    }
  } else {
    if (Lmask) {
      if (Lpicklocal) sector = &SweepLattice2d::sweep_quadrant_mask_picklocal;
      else sector = &SweepLattice2d::sweep_quadrant_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sector = &SweepLattice2d::sweep_quadrant;
    }
  }

  if (Lstrict) ncolor = 4;
  else ncolor = 1;

  // communicator needed between sweep sectors

  comm = new CommLattice2d(spk);
}

/* ---------------------------------------------------------------------- */

SweepLattice2d::~SweepLattice2d()
{
  delete random;
  delete comm;

  memory->destroy_2d_T_array(mask);
  memory->destroy_2d_T_array(ranlat);
}

/* ---------------------------------------------------------------------- */

void SweepLattice2d::init()
{
  applattice = (AppLattice2d *) app;

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

  temperature = applattice->temperature;
  if (temperature != 0.0) t_inverse = 1.0/temperature;
  masklimit = applattice->masklimit;
  
  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;

  nquad = 4;
    
  quad[0].xlo = 1;
  quad[0].xhi = nx_half-1;
  quad[0].ylo = 1;
  quad[0].yhi = ny_half-1;
  
  quad[1].xlo = 1;
  quad[1].xhi = nx_half-1;
  quad[1].ylo = ny_half;
  quad[1].yhi = ny_local;
  
  quad[2].xlo = nx_half;
  quad[2].xhi = nx_local;
  quad[2].ylo = 1;
  quad[2].yhi = ny_half-1;
  
  quad[3].xlo = nx_half;
  quad[3].xhi = nx_local;
  quad[3].ylo = ny_half;
  quad[3].yhi = ny_local;

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    memory->create_2d_T_array(mask,nx_local+2,ny_local+2,
			      "sweeplattice2d:mask");
    
    for (int i = 1; i <= nx_local; i++) 
      for (int j = 1; j <= ny_local; j++) 
	mask[i][j] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    memory->create_2d_T_array(ranlat,nx_local+2,ny_local+2,
			      "sweeplattice2d:ranlat");
    int isite;
    for (int i = 1; i <= nx_local; i++)
      for (int j = 1; j <= ny_local; j++) {
	isite = (i+nx_offset)*ny_global + j + ny_offset;
	ranlat[i][j].init(seed+isite);
      }
  }

  // init communication for ghost sites

  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth);
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice2d::do_sweep(double &dt)
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
   
void SweepLattice2d::sweep_quadrant(int icolor, int iquad)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

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
   
void SweepLattice2d::sweep_quadrant_mask(int icolor, int iquad)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++) mask[i][j] = 0;
  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++) mask[i][j] = 0;

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

      if (efinal <= einitial) continue;
      else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
    }
}

/* ---------------------------------------------------------------------- */

void SweepLattice2d::sweep_quadrant_mask_picklocal(int icolor, int iquad)
{
  int i,j,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++) mask[i][j] = 0;
  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++) mask[i][j] = 0;

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
   
void SweepLattice2d::sweep_quadrant_strict(int icolor, int iquad)
{
  int i,j,i0,j0,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  i0 = (icolor/2 + nx_offset + xlo) % 2;
  j0 = (icolor   + ny_offset + ylo) % 2;

  for (i = xlo+i0; i <= xhi; i += 2)
    for (j = ylo+j0; j <= yhi; j += 2) {
      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);

      newstate = applattice->site_pick_random(i,j,ranlat[i][j].uniform());
      lattice[i][j] = newstate;
      efinal = applattice->site_energy(i,j);

      if (efinal <= einitial) continue;
      else if (temperature == 0.0) lattice[i][j] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse))
	lattice[i][j] = oldstate;
    }
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice2d::sweep_quadrant_mask_strict(int icolor, int iquad)
{
  int i,j,i0,j0,oldstate,newstate;
  double einitial,efinal;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  // unset all masks on boundary
  // may be out of date, due to state change on neighboring processor
  // could reverse comm mask values, but that might be slower

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++) mask[i][j] = 0;
  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++) mask[i][j] = 0;

  i0 = (icolor/2 + nx_offset + xlo) % 2;
  j0 = (icolor   + ny_offset + ylo) % 2;

  // call RNG even if skipping via mask
  // insures same answer no matter where proc boundaries are

  for (i = xlo+i0; i <= xhi; i+=2)
    for (j = ylo+j0; j <= yhi; j+=2) {
      if (mask[i][j]) {
	ranlat[i][j].uniform();
	continue;
      }

      oldstate = lattice[i][j];
      einitial = applattice->site_energy(i,j);
      if (einitial < masklimit) {
	mask[i][j] = 1;
	ranlat[i][j].uniform();
	continue;
      }

      newstate = applattice->site_pick_random(i,j,ranlat[i][j].uniform());
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
