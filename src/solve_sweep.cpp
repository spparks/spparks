/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "solve_sweep.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveSweep::SolveSweep(SPK *spk, int narg, char **arg) :
  Solve(spk,narg,arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
}

/* ---------------------------------------------------------------------- */

SolveSweep::~SolveSweep() {}

/* ---------------------------------------------------------------------- */

void SolveSweep::init(int n, double *propensity) {}

/* ---------------------------------------------------------------------- */

void SolveSweep::update(int n, int *indices, double *propensity) {}

/* ---------------------------------------------------------------------- */

void SolveSweep::update(int n, double *propensity) {}

/* ---------------------------------------------------------------------- */

void SolveSweep::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}

/* ---------------------------------------------------------------------- */

int SolveSweep::event(double *pdt)
{
  return 0;
}
