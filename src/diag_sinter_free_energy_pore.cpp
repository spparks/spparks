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
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_sinter_free_energy_pore.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"
#include "domain.h"

#include <cmath>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagSinterFreeEnergyPore::DiagSinterFreeEnergyPore(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::init()
{
  appsinter = (AppSinter *) app;
  nlocal = appsinter->nlocal;
  interfacialFE = 0.0;
  init_flag = false;
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::compute()
{
	
	if ( !init_flag ) {
		initialize_parameters_calculation();
		init_flag = true;
	}
	
	const int VACANT ( AppSinter::VACANT );
	double total_sites = 0;
	int xgrid, ygrid, zgrid;
	
	int *spin = appsinter->spin;
	tagint *id = appsinter->id;
	int *numneigh = appsinter->numneigh;
	int **neighbor = appsinter->neighbor;
	
	
	double interfacialFEtmp = 0.0;
	for (int i = 0; i < nlocal; i++) {
		appsinter->global_to_grid( id[i], xgrid, ygrid, zgrid );
		if ((xgrid > xstart_) && (xgrid < xend_) &&
			(ygrid > ystart_) && (ygrid < yend_) &&
			(zgrid > zstart_) && (zgrid < zend_) ) {
		  
			total_sites++;
		  
			int ispin = spin[i];
			double surface = 0;
			if ( ispin > VACANT ) { // If I am a grain site add the number of neighbors that are pore sites
				for (int j = 0; j < numneigh[i]; j++)
					if (spin[neighbor[i][j]] == VACANT) surface++;
			}
			interfacialFEtmp += surface;
		}
	}
	
	vector<double> local_info( 2 );
	local_info[0] = interfacialFEtmp;
	local_info[1] = total_sites;
	
	vector<double> info_all( 2, 0 );
	
	MPI_Allreduce(&local_info[0], &info_all[0], 2, MPI_DOUBLE, MPI_SUM, world);
	
	interfacialFE = info_all[0] / info_all[1];	
	
	
//	MPI_Allreduce(&interfacialFEtmp,&interfacialFE,1,MPI_LONG_DOUBLE,MPI_SUM,world);
//  density = appsinter->calculate_density();	
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::stats(char *strtmp)
{
//  sprintf(strtmp," %10.6lf %10g",density,interfacialFE);
  sprintf(strtmp," %10.6lf", interfacialFE);
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::stats_header(char *strtmp)
{
//  sprintf(strtmp," %10s %10s","Density", "InterfFE");
  sprintf(strtmp," %10s", "FE_psite");
}

/* ---------------------------------------------------------------------- */
void DiagSinterFreeEnergyPore::initialize_parameters_calculation()
{
   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;
/*	
	// Determine the central parallelepiped for calculating pore surface
	double gs = appsinter->count_grain_sites();
	double occupied_fraction = (double) gs / (double)( nx * ny * nz );
	double rcube_fraction = pow( occupied_fraction, 1. / 3.);
	int nx_density = (int)floor( nx * rcube_fraction );
	int ny_density = (int)floor( ny * rcube_fraction );
	int nz_density = (int)floor( nz * rcube_fraction );
	
	// Open interval: x > xstart_density and x < xend_density ...
	xstart_ = (nx - nx_density) / 2;
	xend_ = xstart_ + nx_density;
	ystart_ = (ny - ny_density) / 2;
	yend_ = ystart_ + ny_density;
	zstart_ = (nz - nz_density) / 2;
	zend_ = zstart_ + nz_density;
*/	
  int xsize = nx / 3;
  int ysize = ny / 3;
  int zsize = nz / 3;
  
  xstart_ = nx*0.33;
  xend_ = xstart_ + xsize;
  ystart_ = ny*0.33;
  yend_ = ystart_ + ysize;
  zstart_ = nz*0.33;
  zend_ = zstart_ + zsize; 
}

