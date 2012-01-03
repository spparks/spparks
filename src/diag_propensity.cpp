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
#include "diag_propensity.h"
#include "app.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagPropensity::DiagPropensity(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagPropensity::init()
{
  if (!solve)
    error->all(FLERR,"Diag propensity requires KMC solve be performed");

  applattice = (AppLattice *) app;
  nlocal = app->nlocal;
  propensity = 0.0;
}

/* ---------------------------------------------------------------------- */

void DiagPropensity::compute()
{
  applattice->comm->all();

  double ptmp = 0.0;
  for (int i = 0; i < nlocal; i++) ptmp += applattice->site_propensity(i);
  MPI_Allreduce(&ptmp,&propensity,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void DiagPropensity::stats(char *strtmp)
{
  sprintf(strtmp," %10g",propensity);
}

/* ---------------------------------------------------------------------- */

void DiagPropensity::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Propnsty");
}
