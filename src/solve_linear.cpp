/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_linear.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveLinear::SolveLinear(SPK *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  prob = NULL;
}

/* ---------------------------------------------------------------------- */

SolveLinear::~SolveLinear()
{
  delete [] prob;
  delete random;
}

/* ---------------------------------------------------------------------- */

SolveLinear *SolveLinear::clone()
{
  int narg = 2;
  char *arg[2];
  arg[0] = style;
  arg[1] = new char[16];
  sprintf(arg[1],"%d",seed);

  SolveLinear *ptr = new SolveLinear(spk,narg,arg);

  delete [] arg[1];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveLinear::init(int n, double *propensity)
{
  delete [] prob;
  nevents = n;
  prob = new double[n];

  nzeroes = 0;
  sum = 0.0;
  for (int i = 0; i < n; i++) {
    if (propensity[i] == 0.0) nzeroes++;
    prob[i] = propensity[i];
    sum += propensity[i];
  }
}

/* ---------------------------------------------------------------------- */

void SolveLinear::update(int n, int *indices, double *propensity)
{
  int m;
  for (int i = 0; i < n; i++) {
    m = indices[i];
    if (prob[m] == 0.0) nzeroes--;
    if (propensity[m] == 0.0) nzeroes++;
    sum -= prob[m];
    prob[m] = propensity[m];
    sum += propensity[m];
  }
}
/* ---------------------------------------------------------------------- */

void SolveLinear::update(int n, double *propensity)
{
  if (prob[n] == 0.0) nzeroes--;
  if (propensity[n] == 0.0) nzeroes++;
  sum -= prob[n];
  prob[n] = propensity[n];
  sum += propensity[n];
}
/* ---------------------------------------------------------------------- */


void SolveLinear::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}
/* ---------------------------------------------------------------------- */

int SolveLinear::event(double *pdt)
{
  int m;

  if (nzeroes == nevents) {
    sum = 0.0;
    return -1;
  }

  double fraction = sum * random->uniform();
  double partial = 0.0;

  for (m = 0; m < nevents; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
