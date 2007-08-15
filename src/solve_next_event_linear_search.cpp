/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_next_event_linear_search.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveNextEventLinearSearch::
SolveNextEventLinearSearch(SPK *spk, int narg, char **arg) :
  Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
  prob = NULL;
}

/* ---------------------------------------------------------------------- */

SolveNextEventLinearSearch::~SolveNextEventLinearSearch()
{
  delete [] prob;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SolveNextEventLinearSearch::init(int n, double *propensity)
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

void SolveNextEventLinearSearch::update(int n, int *indices, double *propensity)
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

void SolveNextEventLinearSearch::update(int n, double *propensity)
{
  if (prob[n] == 0.0) nzeroes--;
  if (propensity[n] == 0.0) nzeroes++;
  sum -= prob[n];
  prob[n] = propensity[n];
  sum += propensity[n];
}
/* ---------------------------------------------------------------------- */


void SolveNextEventLinearSearch::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}
/* ---------------------------------------------------------------------- */

int SolveNextEventLinearSearch::event(double *pdt)
{
  int m;

  if (nzeroes == nevents) return -1;

  double fraction = sum * random->uniform();
  double partial = 0.0;

  for (m = 0; m < nevents; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
