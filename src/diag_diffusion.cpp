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
#include "string.h"
#include "diag_diffusion.h"
#include "app_diffusion.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagDiffusion::DiagDiffusion(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (strcmp(app->style,"diffusion") != 0)
    error->all(FLERR,"Diag_style diffusion requires app_style diffusion");
}

/* ---------------------------------------------------------------------- */

void DiagDiffusion::init()
{
  appdiff = (AppDiffusion *) app;
  deposit_success = deposit_failed = 0;
}

/* ---------------------------------------------------------------------- */

void DiagDiffusion::compute()
{
  deposit_success = appdiff->ndeposit;
  deposit_failed = appdiff->ndeposit_failed;
  double nfirst = appdiff->nfirst;
  double nsecond = appdiff->nsecond;
  
  MPI_Allreduce(&nfirst,&nfirst_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nsecond,&nsecond_all,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void DiagDiffusion::stats(char *strtmp)
{
  sprintf(strtmp," %10g %10g %10g %10g",
	  deposit_success,deposit_failed,nfirst_all,nsecond_all);
}

/* ---------------------------------------------------------------------- */

void DiagDiffusion::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s %10s %10s",
	  "Deposit","FailedDep","1stHops","2ndHops");
}
