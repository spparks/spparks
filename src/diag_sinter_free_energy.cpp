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
#include "diag_sinter_free_energy.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagSinterFreeEnergy::DiagSinterFreeEnergy(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergy::init()
{
  appsinter = (AppSinter *) app;
  nlocal = appsinter->nlocal;
//  density = 0.0;
  interfacialFE = 0.0;
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergy::compute()
{
//  appsinter->comm->all();
  int *spin = appsinter->spin;
  int *numneigh = appsinter->numneigh;
  int **neighbor = appsinter->neighbor;
  
  const int VACANT ( AppSinter::VACANT );

  double interfacialFEtmp = 0.0;
  for (int i = 0; i < nlocal; i++) {
	int ispin = spin[i];
	double surface = 0;
	if ( ispin > VACANT ) { // If I am a grain site add the number of neighbors that are pore sites
		for (int j = 0; j < numneigh[i]; j++)
			if (spin[neighbor[i][j]] == VACANT) surface++;
	}		
	interfacialFEtmp += surface;
	//interfacialFEtmp += appsinter->site_surface(i);
  }
  MPI_Allreduce(&interfacialFEtmp,&interfacialFE,1,MPI_DOUBLE,MPI_SUM,world);
//  density = appsinter->calculate_density();	
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergy::stats(char *strtmp)
{
//  sprintf(strtmp," %10.6lf %10g",density,interfacialFE);
  sprintf(strtmp," %10Lg", interfacialFE);
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergy::stats_header(char *strtmp)
{
//  sprintf(strtmp," %10s %10s","Density", "InterfFE");
  sprintf(strtmp," %10s", "InterfFE");
}
