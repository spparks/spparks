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
#include "diag_energy2d.h"
#include "app_lattice2d.h"
#include "comm_lattice2d.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEnergy2d::DiagEnergy2d(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style energy2d output file");
	}
      } else {
	error->all("Illegal diag_style energy2d command");
      } 
    } else {
      //      error->all("Illegal diag_style energy2d command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagEnergy2d::~DiagEnergy2d()
{
}

/* ---------------------------------------------------------------------- */

void DiagEnergy2d::init(double time)
{

  if (app->appclass != App::LATTICE2D)
    error->all("Diag style incompatible with app style");

  applattice2d = (AppLattice2d *) app;
  nx_local = applattice2d->nx_local;
  ny_local = applattice2d->ny_local;

  applattice2d->comm->all(applattice2d->lattice);
  
  double etmp = 0.0;
  for (int i = 1; i <= nx_local; i++)
    for (int j = 1; j <= ny_local; j++)
      etmp += applattice2d->site_energy(i,j);
  
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagEnergy2d::compute(double time, int iflag, int done)
{
  double etmp;

  if (diag_delta > 0.0) {
    iflag = check_time(time, done);
  }

  if (iflag) {

    applattice2d->comm->all(applattice2d->lattice);

    etmp = 0.0;
    for (int i = 1; i <= nx_local; i++)
      for (int j = 1; j <= ny_local; j++)
	  etmp += applattice2d->site_energy(i,j);
    
    MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
  }
  
}

/* ---------------------------------------------------------------------- */

void DiagEnergy2d::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy2d::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Energy");
}
