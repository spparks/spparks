/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_grain.h"
#include "sweep_grain.h"
#include "comm_grain_3d.h"
#include "comm_grain_2d.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SweepGrain::SweepGrain(SPK *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  int iarg;

  if (narg < 2) error->all("Illegal sweep_style command");

  // default/dummy values

  appgrain = NULL;

  me = -1;
  nprocs = 0;
  dimension = 0;
  nx_global = 0;
  ny_global = 0;
  nz_global = 0;
  nx_local = 0;
  ny_local = 0;
  nz_local = 0;
  nx_offset = -1;
  ny_offset = -1;
  nz_offset = -1;
  nx_half = -1;
  ny_half = -1;
  nz_half = -1;
  lattice = NULL;
  lat_2d = NULL;
  lat_3d = NULL;
  mask = NULL;
  procwest = -1;
  proceast = -1;
  procsouth = -1;
  procnorth = -1;
  procup = -1;
  procdown = -1;
  nspins = 0;
  temperature = 0.0;
  Lmask = true;
  Lstrict = false;
  nquad = 0;
  masklimit = 0;
  ncolor = 0;

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);

  if (narg == 2) return;
  iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"strict") == 0) {
      iarg++;
      if (iarg == narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg],"yes") == 0) {
	Lstrict = true;
	iarg++;
      } else if (strcmp(arg[iarg],"no") == 0) {
	Lstrict = false;
	iarg++;
      } else error->all("Illegal sweep_style command");
    } else if (strcmp(arg[iarg],"mask") == 0) {
      iarg++;
      if (iarg == narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg],"yes") == 0) {
	Lmask = true;
	iarg++;
      } else if (strcmp(arg[iarg],"no") == 0) {
	Lmask = false;
	iarg++;
      } else error->all("Illegal sweep_style command");
    } else error->all("Illegal sweep_style command");
  }
}

/* ---------------------------------------------------------------------- */

SweepGrain::~SweepGrain()
{
  delete random;
  if (dimension == 2) delete comm_2d;
  else delete comm_3d;

  if (Lmask) {
    memory->destroy_2d_T_array(mask);
  }

  if (Lstrict) {
    memory->destroy_2d_T_array(ranlat);
  }

}

/* ---------------------------------------------------------------------- */

void SweepGrain::init(AppGrain* appgrain_in, const int me_in, const int nprocs_in, const int dimension_in, 
		      const int nx_local_in, const int ny_local_in, const int nz_local_in, 
		      const int nx_global_in, const int ny_global_in, const int nz_global_in, 
		      const int nx_offset_in, const int ny_offset_in, const int nz_offset_in, 
		      int*** lattice_in, 
		      const int procwest_in, const int proceast_in, 
		      const int procsouth_in, const int procnorth_in, 
		      const int procdown_in, const int procup_in, 
		      const int nspins_in, const double temperature_in,
		      int (AppGrain::*es_fp_in)(int,int,int,int,int),
		      void (AppGrain::*um_fp_in)(char***,int,int,int,int))
{
  appgrain = appgrain_in;
  me = me_in;
  nprocs = nprocs_in;
  dimension = dimension_in;
  nx_local = nx_local_in;
  ny_local = ny_local_in;
  nx_global = nx_global_in;
  ny_global = ny_global_in;
  nx_offset = nx_offset_in;
  ny_offset = ny_offset_in;
  lattice = lattice_in;
  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;
  es_fp = es_fp_in;
  um_fp = um_fp_in;

  nx_half = nx_local/2 + 1;
  ny_half = ny_local/2 + 1;
    
  if (dimension == 2) {
    nquad = 4;
    masklimit = 4;
    ncolor = 4;
  } else if (dimension == 3) {
    nquad = 8;
    masklimit = 13;
    ncolor = 8;
    nz_local = nz_local_in;
    nz_global = nz_global_in;
    nz_offset = nz_offset_in;
    procdown = procdown_in;
    procup = procup_in;
    nz_half = nz_local/2 + 1;
  } else {
    error->all("Illegal lattice dimensionality");
  }

  nspins = nspins_in;
  temperature = temperature_in;

  if (dimension == 2) {
    lat_2d = lattice[0];

    comm_2d = new CommGrain2D(spk);
    comm_2d->setup(nx_local,ny_local,nx_half,ny_half,
	      procwest,proceast,procsouth,procnorth);

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
 
    // Initialize mask
    if (Lmask) {
      memory->create_3d_T_array(mask,1,nx_local+2,ny_local+2,
				"sweep_grain:mask");

      for (int i = 0; i <= nx_local+1; i++) 
	for (int j = 0; j <= ny_local+1; j++) 
	  mask[0][i][j] = 0;
    }

    if (Lstrict) {
      int isite;
      // Set up array of random number generators
      memory->create_3d_T_array(ranlat,1,nx_local+2,ny_local+2,
				"sweep_grain:ranlat");
    // construct random number generators
      
      for (int i = 1; i <= nx_local; i++) {
	for (int j = 1; j <= ny_local; j++) {
	  isite = (i+nx_offset)*ny_global + j + ny_offset;
	  ranlat[0][i][j].init(seed+isite);
	}
      }
    }

  } else {
    lat_3d = lattice;

    comm_3d = new CommGrain3D(spk);
    comm_3d->setup(nx_local,ny_local,nz_local,nx_half,ny_half,nz_half,
	      procwest,proceast,procsouth,procnorth,procdown,procup);

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
 
    // Initialize mask
    if (Lmask) {
      memory->create_3d_T_array(mask,nx_local+2,ny_local+2,nz_local+2,
				"sweep_grain:mask");

      for (int i = 0; i <= nx_local+1; i++) 
	for (int j = 0; j <= ny_local+1; j++) 
	  for (int k = 0; k <= nz_local+1; k++) 
	  mask[i][j][k] = 0;
    }

    if (Lstrict) {
      int isite;
      // Set up array of random number generators
      memory->create_3d_T_array(ranlat,nx_local+2,ny_local+2,nz_local+2,
				"sweep_grain:ranlat");
    // construct random number generators
      
      for (int i = 1; i <= nx_local; i++) {
	for (int j = 1; j <= ny_local; j++) {
	  for (int k = 1; k <= nz_local; k++) {
	    isite = (i+nx_offset)*ny_global*nz_global + 
	      (j+ny_offset)*nz_global + k + nz_offset;
	    ranlat[i][j][k].init(seed+isite);
	  }
	}
      }
    }

  }

}

/* ----------------------------------------------------------------------
   perform one sweep over domain
 ------------------------------------------------------------------------- */

void SweepGrain::do_sweep()
{
  if (Lstrict) {
    // Loop over the lattice colors a la checkerboarding for 2D Ising model
    for (int icolor = 0; icolor < ncolor; icolor++) {
      // Loop over the sectors
      for (int iquad = 0; iquad < nquad; iquad++) {
	timer->stamp();
	if (dimension == 2) comm_2d->communicate(lat_2d,iquad);
	else comm_3d->communicate(lat_3d,iquad);
	timer->stamp(TIME_COMM);

	timer->stamp();
	if (Lmask) {
	  // Eventually hide these variations using functions pointers to AppGrain methods
	  if (dimension == 2) sweep_quadrant_mask_strict_2d(icolor, iquad);
	  else sweep_quadrant_mask_strict_3d(icolor, iquad);
	} else {
	  if (dimension == 2) sweep_quadrant_strict_2d(icolor, iquad);
	  else sweep_quadrant_strict_3d(icolor, iquad);
	}
	timer->stamp(TIME_SOLVE);
      }
    }
  } else {
    // Loop over the sectors
    for (int iquad = 0; iquad < nquad; iquad++) {
      timer->stamp();
      if (dimension == 2) comm_2d->communicate(lat_2d,iquad);
      else comm_3d->communicate(lat_3d,iquad);
      timer->stamp(TIME_COMM);
	
      timer->stamp();
      if (Lmask) {
	if (dimension == 2) sweep_quadrant_mask_2d(iquad);
	else sweep_quadrant_mask_3d(iquad);
      } else {
	if (dimension == 2) sweep_quadrant_2d(iquad);
	else sweep_quadrant_3d(iquad);
      }
      timer->stamp(TIME_SOLVE);
    }
  }
    
}
    

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_2d(int iquad)
{
  int i,j,iold,inew,nold,nnew;
  
  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  for (i = xlo; i <= xhi; i++) {
    for (j = ylo; j <= yhi; j++) {
      nold = (appgrain->*es_fp)(lat_2d[i][j],i,j,0,0);

      inew = random->irandom(nspins);
      nnew = (appgrain->*es_fp)(inew,i,j,0,0);

      if (nold >= nnew) {
	lat_2d[i][j] = inew;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_mask_2d(int iquad)
{
  int i,j,iold,inew,nold,nnew;
  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  // Unset all masks at boundary, as they may be out of date,
  // due to spin change on neighboring processor.
  // Alternatively, we could reverse communicate the masks, but that
  // might work out to be even slower....

  char** mask_2d_tmp = mask[0];
  char*** mask_3d_tmp = mask;

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++)
    mask_2d_tmp[i][j] = 0;
  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++)
    mask_2d_tmp[i][j] = 0;

  for (i = xlo; i <= xhi; i++) {
    for (j = ylo; j <= yhi; j++) {
      // Check if mask is set 
      if (mask_2d_tmp[i][j]) {
	continue;
      }

      nold = (appgrain->*es_fp)(lat_2d[i][j],i,j,0,0);
      // Check if mask can be set
      if (nold < masklimit) {
	mask_2d_tmp[i][j] = 1;
	continue;
      }

      inew = random->irandom(nspins);
      nnew = (appgrain->*es_fp)(inew,i,j,0,0);

      if (nold >= nnew) {
	lat_2d[i][j] = inew;
	(appgrain->*um_fp)(mask_3d_tmp,i,j,0,0);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_strict_2d(int icolor, int iquad)
{
  int i,j,iold,inew,nold,nnew,i0,j0;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  RandomPark** ranlat_2d_tmp = ranlat[0];
  i0 = (icolor/2 + nx_offset + xlo) % 2;
  j0 = (icolor   + ny_offset + ylo) % 2;

//   printf("icolor = %d iquad = %d \n",icolor,iquad);
//   printf("xlo = %d i0 = %d xhi = %d \n",xlo,i0,xhi);
//   printf("ylo = %d j0 = %d yhi = %d \n",ylo,j0,yhi);

  for (i = xlo+i0; i <= xhi; i+=2) {
    for (j = ylo+j0; j <= yhi; j+=2) {
      
      nold = (appgrain->*es_fp)(lat_2d[i][j],i,j,0,0);
      
      inew = ranlat_2d_tmp[i][j].irandom(nspins);
      nnew = (appgrain->*es_fp)(inew,i,j,0,0);
      
      if (nold >= nnew) {
	lat_2d[i][j] = inew;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_mask_strict_2d(int icolor, int iquad)
{
  int i,j,iold,inew,nold,nnew,i0,j0;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;

  // Unset all masks at boundary, as they may be out of date,
  // due to spin change on neighboring processor.
  // Alternatively, we could reverse communicate the masks, but that
  // might work out to be even slower....

  char** mask_2d_tmp = mask[0];
  char*** mask_3d_tmp = mask;

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++)
    mask_2d_tmp[i][j] = 0;
  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++)
    mask_2d_tmp[i][j] = 0;

  RandomPark** ranlat_2d_tmp = ranlat[0];
  i0 = (icolor/2 + nx_offset + xlo) % 2;
  j0 = (icolor   + ny_offset + ylo) % 2;

  for (i = xlo+i0; i <= xhi; i+=2) {
    for (j = ylo+j0; j <= yhi; j+=2) {
      // Check if mask is set 
      if (mask_2d_tmp[i][j]) {
	ranlat_2d_tmp[i][j].irandom(nspins);
	continue;
      }
      
      nold = (appgrain->*es_fp)(lat_2d[i][j],i,j,0,0);
      // Check if mask can be set
      if (nold < masklimit) {
	mask_2d_tmp[i][j] = 1;
	ranlat_2d_tmp[i][j].irandom(nspins);
	continue;
      }
      
      inew = ranlat_2d_tmp[i][j].irandom(nspins);
      nnew = (appgrain->*es_fp)(inew,i,j,0,0);
      
      if (nold >= nnew) {
	lat_2d[i][j] = inew;
	(appgrain->*um_fp)(mask_3d_tmp,i,j,0,0);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_3d(int iquad)
{
  int i,j,k,iold,inew,nold,nnew;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  for (i = xlo; i <= xhi; i++) {
    for (j = ylo; j <= yhi; j++) {
      for (k = zlo;  k <= zhi; k++) {
	nold = (appgrain->*es_fp)(lat_3d[i][j][k],i,j,k,0);

	inew = random->irandom(nspins);
	nnew = (appgrain->*es_fp)(inew,i,j,k,0);

	if (nold >= nnew) {
	  lat_3d[i][j][k] = inew;
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_mask_3d(int iquad)
{
  int i,j,k,iold,inew,nold,nnew;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  // Unset all masks at boundary, as they may be out of date,
  // due to spin change on neighboring processor.
  // Alternatively, we could reverse communicate the masks, but that
  // might work out to be even slower....

  char*** mask_3d_tmp = mask;
  if (zlo == 1) k = zlo;
  else if (zhi == nz_local) k = zhi;
  else error->one("Failed to find z-boundary in sweep_quadrant_maks_3d()");
  for (i = xlo; i <= xhi; i++)
  for (j = ylo; j <= yhi; j++)
    mask_3d_tmp[i][j][k] = 0;

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++)
  for (k = zlo; k <= zhi; k++)
    mask_3d_tmp[i][j][k] = 0;

  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++)
  for (k = zlo; k <= zhi; k++)
    mask_3d_tmp[i][j][k] = 0;

  for (i = xlo; i <= xhi; i++) {
    for (j = ylo; j <= yhi; j++) {
      for (k = zlo;  k <= zhi; k++) {
	// Check if mask is set 
	if (mask_3d_tmp[i][j][k]) {
	  continue;
	}

	nold = (appgrain->*es_fp)(lat_3d[i][j][k],i,j,k,0);
	// Check if mask can be set
	if (nold < masklimit) {
	  mask_3d_tmp[i][j][k] = 1;
	  continue;
	}

	inew = random->irandom(nspins);
	nnew = (appgrain->*es_fp)(inew,i,j,k,0);

	if (nold >= nnew) {
	  lat_3d[i][j][k] = inew;
	  (appgrain->*um_fp)(mask_3d_tmp,i,j,k,0);
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_strict_3d(int icolor, int iquad)
{
  int i,j,k,iold,inew,nold,nnew,i0,j0,k0;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  i0 = (icolor/4 + nx_offset + xlo) % 2;
  j0 = (icolor/2 + ny_offset + ylo) % 2;
  k0 = (icolor   + nz_offset + zlo) % 2;

//   printf("icolor = %d iquad = %d \n",icolor,iquad);
//   printf("xlo = %d i0 = %d xhi = %d \n",xlo,i0,xhi);
//   printf("ylo = %d j0 = %d yhi = %d \n",ylo,j0,yhi);

  for (i = xlo+i0; i <= xhi; i+=2) {
    for (j = ylo+j0; j <= yhi; j+=2) {
      for (k = zlo+k0;  k <= zhi; k+=2) {
	nold = (appgrain->*es_fp)(lat_3d[i][j][k],i,j,k,0);
      
	inew = ranlat[i][j][k].irandom(nspins);
	nnew = (appgrain->*es_fp)(inew,i,j,k,0);
	
	if (nold >= nnew) {
	  lat_3d[i][j][k] = inew;
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
   
void SweepGrain::sweep_quadrant_mask_strict_3d(int icolor, int iquad)
{
  int i,j,k,iold,inew,nold,nnew,i0,j0,k0;

  int xlo = quad[iquad].xlo;
  int xhi = quad[iquad].xhi;
  int ylo = quad[iquad].ylo;
  int yhi = quad[iquad].yhi;
  int zlo = quad[iquad].zlo;
  int zhi = quad[iquad].zhi;

  // Unset all masks at boundary, as they may be out of date,
  // due to spin change on neighboring processor.
  // Alternatively, we could reverse communicate the masks, but that
  // might work out to be even slower....

  char*** mask_3d_tmp = mask;
  if (zlo == 1) k = zlo;
  else if (zhi == nz_local) k = zhi;
  else error->one("Failed to find z-boundary in sweep_quadrant_maks_3d()");
  for (i = xlo; i <= xhi; i++)
  for (j = ylo; j <= yhi; j++)
    mask_3d_tmp[i][j][k] = 0;

  if (ylo == 1) j = ylo;
  else if (yhi == ny_local) j = yhi;
  for (i = xlo; i <= xhi; i++)
  for (k = zlo; k <= zhi; k++)
    mask_3d_tmp[i][j][k] = 0;

  if (xlo == 1) i = xlo;
  else if (xhi == nx_local) i = xhi;
  for (j = ylo; j <= yhi; j++)
  for (k = zlo; k <= zhi; k++)
    mask_3d_tmp[i][j][k] = 0;

  i0 = (icolor/4 + nx_offset + xlo) % 2;
  j0 = (icolor/2 + ny_offset + ylo) % 2;
  k0 = (icolor   + nz_offset + zlo) % 2;

  for (i = xlo+i0; i <= xhi; i+=2) {
    for (j = ylo+j0; j <= yhi; j+=2) {
      for (k = zlo+k0;  k <= zhi; k+=2) {
	// Check if mask is set 
	if (mask_3d_tmp[i][j][k]) {
	  ranlat[i][j][k].irandom(nspins);
	  continue;
	}
      
	nold = (appgrain->*es_fp)(lat_3d[i][j][k],i,j,k,0);
	// Check if mask can be set
	if (nold < masklimit) {
	  mask_3d_tmp[i][j][k] = 1;
	  ranlat[i][j][k].irandom(nspins);
	  continue;
	}
      
	inew = ranlat[i][j][k].irandom(nspins);
	nnew = (appgrain->*es_fp)(inew,i,j,k,0);
	
	if (nold >= nnew) {
	  lat_3d[i][j][k] = inew;
	  (appgrain->*um_fp)(mask_3d_tmp,i,j,k,0);
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Compute total energy of system, quadrant by quandrant
------------------------------------------------------------------------- */

double SweepGrain::compute_energy()
{
  double energy_local,energy_global,energy_quad;
  energy_local = 0.0;

  for (int iquad = 0; iquad < nquad; iquad++) {

    if (dimension == 2) comm_2d->communicate(lat_2d,iquad);
    else comm_3d->communicate(lat_3d,iquad);

    energy_quad = energy_quadrant(iquad);
    energy_local += energy_quad;
  }

  MPI_Allreduce(&energy_local,&energy_global,1,MPI_DOUBLE,MPI_SUM,world);
  
  return energy_global;
}

/* ----------------------------------------------------------------------
   Compute energy in one quadrant
------------------------------------------------------------------------- */

double SweepGrain::energy_quadrant(const int iquad)
{
  int i,j,k;
  double energy;

  energy = 0.0;


  if (dimension == 2) {
    for (i = quad[iquad].xlo; i <= quad[iquad].xhi; i++) {
      for (j = quad[iquad].ylo; j <= quad[iquad].yhi; j++) {
	energy+=(appgrain->*es_fp)(lat_2d[i][j],i,j,0,0);
      }
    }
  } else {
    for (i = quad[iquad].xlo; i <= quad[iquad].xhi; i++) {
      for (j = quad[iquad].ylo; j <= quad[iquad].yhi; j++) {
	for (k = quad[iquad].zlo; k <= quad[iquad].zhi; k++) {
	  energy+=(appgrain->*es_fp)(lat_3d[i][j][k],i,j,k,0);
	}
      }
    }
  }

  return energy;

}
