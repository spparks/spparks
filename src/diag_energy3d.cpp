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
#include "diag_energy3d.h"
#include "app_lattice3d.h"
#include "comm_lattice3d.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEnergy3d::DiagEnergy3d(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style energy3d output file");
	}
      } else {
	error->all("Illegal diag_style energy3d command");
      } 
    } else {
      //      error->all("Illegal diag_style energy3d command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagEnergy3d::~DiagEnergy3d()
{
}

/* ---------------------------------------------------------------------- */

void DiagEnergy3d::init(double time)
{
  if (app->appclass != App::LATTICE3D)
    error->all("Diag style incompatible with app style");

  applattice3d = (AppLattice3d *) app;
  nx_local = applattice3d->nx_local;
  ny_local = applattice3d->ny_local;
  nz_local = applattice3d->nz_local;

  applattice3d->comm->all(applattice3d->lattice);
  
  double etmp = 0.0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++)
      for (int k = 1; k <= nz_local; k++)
	etmp += applattice3d->site_energy(i,j,k);
  
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagEnergy3d::compute(double time, int iflag, int done)
{
  double etmp;

  if (diag_delta > 0.0) {
    iflag = check_time(time, done);
  }

  if (iflag || done) {

    applattice3d->comm->all(applattice3d->lattice);

    etmp = 0.0;
    for (int i = 1; i <= nx_local; i++)
      for (int j = 1; j <= ny_local; j++)
	for (int k = 1; k <= nz_local; k++)
	  etmp += applattice3d->site_energy(i,j,k);
    
    MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
  }
  
}

/* ---------------------------------------------------------------------- */

void DiagEnergy3d::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy3d::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Energy");
}
