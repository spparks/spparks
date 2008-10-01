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
#include "diag_energy.h"
#include "app_lattice.h"
#include "comm_lattice.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEnergy::DiagEnergy(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style energy output file");
	}
      } else {
	error->all("Illegal diag_style energy command");
      } 
    } else {
      //      error->all("Illegal diag_style energy command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagEnergy::~DiagEnergy()
{
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::init(double time)
{
  if (app->appclass != App::LATTICE)
    error->all("Diag style incompatible with app style");

  applattice = (AppLattice *) app;

  nlocal = applattice->nlocal;

  applattice->comm->all();
  
  double etmp = 0.0;
  for (int i = 0; i < nlocal; i++)
    etmp += applattice->site_energy(i);
  
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagEnergy::compute(double time, int iflag, int done)
{
  double etmp;

  if (stats_flag == 0) {
    iflag = check_time(time, done);
  }

  if (iflag) {
    applattice->comm->all();

    etmp = 0.0;
    for (int i = 0; i < nlocal; i++) 
      etmp += applattice->site_energy(i);

    
    MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  }
  
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Energy");
}
