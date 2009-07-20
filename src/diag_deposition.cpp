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
#include "diag_deposition.h"
#include "app_diffusion2.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagDeposition::DiagDeposition(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (strcmp(app->style,"diffusion2") != 0)
    error->all("Diag_style deposition requires app_style diffusion");
}

/* ---------------------------------------------------------------------- */

void DiagDeposition::init()
{
  appdiff = (AppDiffusion2 *) app;
  deposit_success = deposit_failed = 0;
}

/* ---------------------------------------------------------------------- */

void DiagDeposition::compute()
{
  deposit_success = appdiff->ndeposit;
  deposit_failed = appdiff->ndeposit_failed;
}

/* ---------------------------------------------------------------------- */

void DiagDeposition::stats(char *strtmp)
{
  sprintf(strtmp,"%10g %10g",deposit_success,deposit_failed);
}

/* ---------------------------------------------------------------------- */

void DiagDeposition::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s","Deposit","FailedDep");
}
