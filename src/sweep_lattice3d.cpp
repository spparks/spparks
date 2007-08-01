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

#define SolveInclude
#include "style.h"
#undef SolveInclude

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SweepLattice3d::SweepLattice3d(SPK *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  int seed = atoi(arg[1]);
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
    sector = &SweepLattice3d::sweep_quadrant_kmc;
  } else if (Lstrict) {
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

  nquad = 8;
  for (int i = 0; i < nquad; i++) {
    quad[i].propensity = NULL;
    quad[i].site2ijk = NULL;
    quad[i].sites = NULL;
  }

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

  if (Lkmc) {
    for (int iquad = 0; iquad < nquad; iquad++) {
      delete quad[iquad].solve;
      memory->sfree(quad[iquad].propensity);
      memory->destroy_2d_T_array(quad[iquad].site2ijk);
      memory->sfree(quad[iquad].sites);
    }
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
  masklimit = applattice->masklimit;
  
  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

  quad[0].xlo = 1;
  quad[0].xhi = nx_half-1;
  quad[0].ylo = 1;
  quad[0].yhi = ny_half-1;
  quad[0].zlo = 1;
  quad[0].zhi = nz_half-1;
  quad[0].nx = quad[0].xhi - quad[0].xlo + 1;
  quad[0].ny = quad[0].yhi - quad[0].ylo + 1;
  quad[0].nz = quad[0].zhi - quad[0].zlo + 1;
  memory->sfree(quad[0].propensity);
  memory->destroy_2d_T_array(quad[0].site2ijk);
  memory->sfree(quad[0].sites);
  
  quad[1].xlo = 1;
  quad[1].xhi = nx_half-1;
  quad[1].ylo = 1;
  quad[1].yhi = ny_half-1;
  quad[1].zlo = nz_half;
  quad[1].zhi = nz_local;
  quad[1].nx = quad[1].xhi - quad[1].xlo + 1;
  quad[1].ny = quad[1].yhi - quad[1].ylo + 1;
  quad[1].nz = quad[1].zhi - quad[1].zlo + 1;
  memory->sfree(quad[1].propensity);
  memory->destroy_2d_T_array(quad[1].site2ijk);
  memory->sfree(quad[1].sites);

  quad[2].xlo = 1;
  quad[2].xhi = nx_half-1;
  quad[2].ylo = ny_half;
  quad[2].yhi = ny_local;
  quad[2].zlo = 1;
  quad[2].zhi = nz_half-1;
  quad[2].nx = quad[2].xhi - quad[2].xlo + 1;
  quad[2].ny = quad[2].yhi - quad[2].ylo + 1;
  quad[2].nz = quad[2].zhi - quad[2].zlo + 1;
  memory->sfree(quad[2].propensity);
  memory->destroy_2d_T_array(quad[2].site2ijk);
  memory->sfree(quad[2].sites);

  quad[3].xlo = 1;
  quad[3].xhi = nx_half-1;
  quad[3].ylo = ny_half;
  quad[3].yhi = ny_local;
  quad[3].zlo = nz_half;
  quad[3].zhi = nz_local;
  quad[3].nx = quad[3].xhi - quad[3].xlo + 1;
  quad[3].ny = quad[3].yhi - quad[3].ylo + 1;
  quad[3].nz = quad[3].zhi - quad[3].zlo + 1;
  memory->sfree(quad[3].propensity);
  memory->destroy_2d_T_array(quad[3].site2ijk);
  memory->sfree(quad[3].sites);

  quad[4].xlo = nx_half;
  quad[4].xhi = nx_local;
  quad[4].ylo = 1;
  quad[4].yhi = ny_half-1;
  quad[4].zlo = 1;
  quad[4].zhi = nz_half-1;
  quad[4].nx = quad[4].xhi - quad[4].xlo + 1;
  quad[4].ny = quad[4].yhi - quad[4].ylo + 1;
  quad[4].nz = quad[4].zhi - quad[4].zlo + 1;
  memory->sfree(quad[4].propensity);
  memory->destroy_2d_T_array(quad[4].site2ijk);
  memory->sfree(quad[4].sites);

  quad[5].xlo = nx_half;
  quad[5].xhi = nx_local;
  quad[5].ylo = 1;
  quad[5].yhi = ny_half-1;
  quad[5].zlo = nz_half;
  quad[5].zhi = nz_local;
  quad[5].nx = quad[5].xhi - quad[5].xlo + 1;
  quad[5].ny = quad[5].yhi - quad[5].ylo + 1;
  quad[5].nz = quad[5].zhi - quad[5].zlo + 1;
  memory->sfree(quad[5].propensity);
  memory->destroy_2d_T_array(quad[5].site2ijk);
  memory->sfree(quad[5].sites);

  quad[6].xlo = nx_half;
  quad[6].xhi = nx_local;
  quad[6].ylo = ny_half;
  quad[6].yhi = ny_local;
  quad[6].zlo = 1;
  quad[6].zhi = nz_half-1;
  quad[6].nx = quad[6].xhi - quad[6].xlo + 1;
  quad[6].ny = quad[6].yhi - quad[6].ylo + 1;
  quad[6].nz = quad[6].zhi - quad[6].zlo + 1;
  memory->sfree(quad[6].propensity);
  memory->destroy_2d_T_array(quad[6].site2ijk);
  memory->sfree(quad[6].sites);

  quad[7].xlo = nx_half;
  quad[7].xhi = nx_local;
  quad[7].ylo = ny_half;
  quad[7].yhi = ny_local;
  quad[7].zlo = nz_half;
  quad[7].zhi = nz_local;
  quad[7].nx = quad[7].xhi - quad[7].xlo + 1;
  quad[7].ny = quad[7].yhi - quad[7].ylo + 1;
  quad[7].nz = quad[7].zhi - quad[7].zlo + 1;
  memory->sfree(quad[7].propensity);
  memory->destroy_2d_T_array(quad[7].site2ijk);
  memory->sfree(quad[7].sites);

  // init communication for ghost sites

  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup);

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

  // setup KMC solver and propensity arrays for each quadrant
  // set ijk2site and site2ijk values in app to reflect quadrant mapping
  // propensity init requires ghost cell info for entire sub-domain

  if (Lkmc) {
    int i,j,k,m;

    comm->all(lattice);

    char *arg[2];
    arg[0] = solve->style;
    arg[1] = "12345";           // this line is a kludge

    for (int iquad = 0; iquad < nquad; iquad++) {
      if (strcmp(arg[0],"none") == 0) solve = NULL;

#define SolveClass
#define SolveStyle(key,Class) \
      else if (strcmp(arg[0],#key) == 0) quad[iquad].solve = new Class(spk,2,arg);
#include "style.h"
#undef SolveClass

      int nsites = quad[iquad].nx * quad[iquad].ny * quad[iquad].nz;
      int nborder = 2 * (quad[iquad].nx*quad[iquad].ny + 
			 quad[iquad].ny*quad[iquad].nz + 
			 quad[iquad].nx*quad[iquad].nz );
      quad[iquad].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep:propensity");
      memory->create_2d_T_array(quad[iquad].site2ijk,nsites,3,
				"sweep:ijk2site");
      quad[iquad].sites =
	(int*) memory->smalloc(nborder*sizeof(int),"sweep:sites");

      for (i = quad[iquad].xlo; i <= quad[iquad].xhi; i++)
	for (j = quad[iquad].ylo; j <= quad[iquad].yhi; j++)
	  for (k = quad[iquad].zlo; k <= quad[iquad].zhi; k++)
	    ijk2site[i][j][k] = 
	      (i-quad[iquad].xlo)*quad[iquad].ny*quad[iquad].nz + 
	      (j-quad[iquad].ylo)*quad[iquad].nz + k-quad[iquad].zlo;

      for (m = 0; m < nsites; m++) {
	i = m / quad[iquad].ny/quad[iquad].nz + 1;
	j = (m / quad[iquad].nz) % quad[iquad].ny + 1;
	k = m % quad[iquad].nz + 1;
	quad[iquad].site2ijk[m][0] = i + quad[iquad].xlo - 1;
	quad[iquad].site2ijk[m][1] = j + quad[iquad].ylo - 1;
	quad[iquad].site2ijk[m][2] = k + quad[iquad].zlo - 1;
      }

      for (i = quad[iquad].xlo; i <= quad[iquad].xhi; i++)
	for (j = quad[iquad].ylo; j <= quad[iquad].yhi; j++)
	  for (k = quad[iquad].zlo; k <= quad[iquad].zhi; k++)
	    quad[iquad].propensity[ijk2site[i][j][k]] = 
	      applattice->site_propensity(i,j,k,0);

      quad[iquad].solve->init(nsites,quad[iquad].propensity);
    }
  }
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

/* ----------------------------------------------------------------------
   generate events in each quadrant using KMC solver
   sweep over one sector of sites
   skip sites that can't change via mask
 ------------------------------------------------------------------------- */

void SweepLattice3d::sweep_quadrant_kmc(int icolor, int iquad)
{
  double dt,time;
  int done,isite,i,j,k;

  // extract sector specific info from quad struct

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  Solve *solve = quad[iquad].solve;
  double *propensity = quad[iquad].propensity;
  int **site2ijk = quad[iquad].site2ijk;
  int *sites = quad[iquad].sites;

  // temporarily reset values in applattice
  // sector bounds, propensity array, solver

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

  // update propensities on all 6 sector faces
  // necessary since ghosts of sector may have changed

  int nsites = 0;

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

  solve->update(nsites,sites,propensity);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
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
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = hold_propensity;
}
