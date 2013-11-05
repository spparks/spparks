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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_linear.h"
#include "domain.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

SolveLinear::SolveLinear(SPPARKS *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 1) error->all(FLERR,"Illegal solve command");

  // each proc uses different initial RNG seed

  random = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  random->reset(seed,spk->domain->me,100);

  prob = NULL;
}

/* ---------------------------------------------------------------------- */

SolveLinear::~SolveLinear()
{
  delete random;
  memory->destroy(prob);
}

/* ---------------------------------------------------------------------- */

SolveLinear *SolveLinear::clone()
{
  int narg = 1;
  char *arg[1];
  arg[0] = style;

  SolveLinear *ptr = new SolveLinear(spk,narg,arg);

  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveLinear::init(int n, double *propensity)
{
  delete [] prob;
  nevents = n;
  memory->create(prob,n,"solve/linear:prob");

  sum = 0.0;
  num_active = 0;

  for (int i = 0; i < n; i++) {
    if (propensity[i] > 0.0) num_active++;
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
    if (prob[m] > 0.0) num_active--;
    if (propensity[m] > 0.0) num_active++;
    sum -= prob[m];
    prob[m] = propensity[m];
    sum += propensity[m];
  }
}
/* ---------------------------------------------------------------------- */

void SolveLinear::update(int n, double *propensity)
{
  if (prob[n] > 0.0) num_active--;
  if (propensity[n] > 0.0) num_active++;
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

  if (num_active == 0) {
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

  if (m < nevents) return m;
  return nevents-1;
}

