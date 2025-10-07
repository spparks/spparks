/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
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
