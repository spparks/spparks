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
#include "diag_energy.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "comm_lattice.h"
#include "comm_off_lattice.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagEnergy::DiagEnergy(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass == App::LATTICE) latticeflag = 1;
  else if (app->appclass == App::OFF_LATTICE) latticeflag = 0;
  else error->all(FLERR,"Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::init()
{
  if (latticeflag) applattice = (AppLattice *) app;
  else appofflattice = (AppOffLattice *) app;

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::compute()
{
  int nlocal = app->nlocal;
  if (latticeflag) applattice->comm->all();
  else appofflattice->comm->all();

  double etmp = 0.0;
  if (latticeflag)
    for (int i = 0; i < nlocal; i++) etmp += applattice->site_energy(i);
  else
    for (int i = 0; i < nlocal; i++) etmp += appofflattice->site_energy(i);

  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats(char *strtmp)
{
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Energy");
}
