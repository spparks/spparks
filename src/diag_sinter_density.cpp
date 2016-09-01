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
#include "diag_sinter_density.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"
#include "domain.h"

#include <cmath>

using namespace SPPARKS_NS;


/* ---------------------------------------------------------------------- */

DiagSinterDensity::DiagSinterDensity(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::init()
{
  appsinter = (AppSinter *) app;
  nlocal = appsinter->nlocal;
  density = 0.0; 

}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::compute()
{
//  applattice->comm->all();

	double grain_sites = 0;
	double total_sites = 0;
	int xgrid, ygrid, zgrid;
	static bool init_flag = false;
	
	if ( !init_flag ) {
		initialize_parameters_density_calculation();
		init_flag = true;
//		printf("xstart_den: %d xend_den: %d ystart_den: %d yend_den: %d zstart_den: %d zend_den: %d\n", xstart_density, xend_density, ystart_density, yend_density, zstart_density, zend_density ); 
	}
	
	int *spin = appsinter->spin;
	tagint *id = appsinter->id;
	const int VACANT ( AppSinter::VACANT );
	
	for ( int i = 0; i < nlocal; i++ ) {
		appsinter->global_to_grid( id[i], xgrid, ygrid, zgrid );
		if (	(xgrid > xstart_density) && (xgrid < xend_density) &&
				(ygrid > ystart_density) && (ygrid < yend_density) &&
				(zgrid > zstart_density) && (zgrid < zend_density) ) {
			total_sites++;
			if ( spin[i] > VACANT ) {
				grain_sites++;
			}
		}
	}
	vector<double> local_density_info( 2 );
	local_density_info[0] = grain_sites;
	local_density_info[1] = total_sites;
	
//	printf( "grain sites: %d total sites: %d\n", (int)grain_sites, (int)total_sites );
	
	vector<double> density_info_all( 2, 0 );
	
	MPI_Allreduce(&local_density_info[0], &density_info_all[0], 2, MPI_DOUBLE, MPI_SUM, world);
	
	density = density_info_all[0] / density_info_all[1];


//  density = appsinter->calculate_density();	
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::stats(char *strtmp)
{
  sprintf(strtmp," %10.6lf",density);
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Density");
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::initialize_parameters_density_calculation()
{

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;
/*
	// Determine the central parallelepiped for calculating density
	double gs = appsinter->count_grain_sites();
	double occupied_fraction = (double) gs / (double)( nx * ny * nz );
	double rcube_fraction = pow( occupied_fraction, 1. / 3.);
	int nx_density = (int)floor( nx * rcube_fraction );
	int ny_density = (int)floor( ny * rcube_fraction );
	int nz_density = (int)floor( nz * rcube_fraction );
  
	// Open interval: x > xstart_density and x < xend_density ...
	xstart_density = (nx - nx_density) / 2;
	xend_density = xstart_density + nx_density;
	ystart_density = (ny - ny_density) / 2;
	yend_density = ystart_density + ny_density;
	zstart_density = (nz - nz_density) / 2;
	zend_density = zstart_density + nz_density;
*/
  int xsize = nx / 3;
  int ysize = ny / 3;
  int zsize = nz / 3;
  
  xstart_density = nx*0.33;
  xend_density = xstart_density + xsize;
  ystart_density = ny*0.33;
  yend_density = ystart_density + ysize;
  zstart_density = nz*0.33;
  zend_density = zstart_density + zsize; 
}
