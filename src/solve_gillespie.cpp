/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_gillespie.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveGillespie::SolveGillespie(SPK *spk, int narg, char **arg) : 
  Solve(spk,narg,arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
  prob = NULL;
}

/* ---------------------------------------------------------------------- */

SolveGillespie::~SolveGillespie()
{
  delete random;
  delete [] prob;
}

/* ---------------------------------------------------------------------- */

void SolveGillespie::init(int n, double *propensity)
{
  delete [] prob;
  nreactions = n;
  prob = new double[n];

  sum = 0.0;
  for (int i = 0; i < n; i++) {
    prob[i] = propensity[i];
    sum += propensity[i];
  }
}

/* ---------------------------------------------------------------------- */

void SolveGillespie::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    sum -= prob[indices[i]];
    prob[indices[i]] = propensity[i];
    sum += propensity[i];
  }
}

/* ---------------------------------------------------------------------- */

int SolveGillespie::event(double *pdt)
{
  int m;

  if (sum == 0.0) return -1;
  double fraction = sum * random->uniform();
  
  double partial = 0.0;
  for (m = 0; m < nreactions; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
