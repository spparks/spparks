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
#include "solve_group.h"
#include "groups.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SolveGroup::SolveGroup(SPPARKS *spk, int narg, char **arg) :
  Solve(spk, narg, arg)
{
  if (narg < 3) error->all("Illegal solve command");
  
  hi = atof(arg[1]);
  lo = atof(arg[2]);

  if (lo <= 0.0 || lo >= hi)
    error->all("Invalid probability bounds for solve_style group");

  ngroups = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ngroup") == 0) {
      if (iarg+2 > narg) error->all("Illegal solve_style group command");
      ngroups = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal solve_style group command");
  }

  random = new RandomPark(ranmaster->uniform());
  groups = new Groups(spk,hi,lo,ngroups);
  p = NULL;
}

/* ---------------------------------------------------------------------- */

SolveGroup::~SolveGroup()
{
  delete random;
  delete groups;
  delete [] p;
}

/* ---------------------------------------------------------------------- */

SolveGroup *SolveGroup::clone()
{
  int narg = 5;
  char *arg[5];

  arg[0] = style;
  arg[1] = new char[16];
  arg[2] = new char[16];
  arg[3] = "ngroup";
  arg[4] = new char[16];

  sprintf(arg[1],"%g",hi);
  sprintf(arg[2],"%g",lo);
  sprintf(arg[4],"%d",ngroups);

  SolveGroup *ptr = new SolveGroup(spk,narg,arg);

  delete [] arg[1];
  delete [] arg[2];
  delete [] arg[4];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveGroup::init(int n, double *propensity)
{
  nevents = n;
  num_active = 0;

  delete [] p;
  p = new double[n];

  sum = 0.0;
  for (int i = 0; i < n; i++) {
    if (propensity[i] > 0.0) num_active++;
    p[i] = MAX(propensity[i],lo);
    p[i] = MIN(p[i],hi);
    sum += p[i];
  }

  groups->partition_init(p,n);
}

/* ---------------------------------------------------------------------- */

void SolveGroup::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    int j = indices[i];
    double pt = propensity[j];
    if (p[j] != pt) {
      if (p[j] > 0.0) num_active--;
      if (pt > 0.0) num_active++;
      sum -= p[j];
      groups->alter_element(j,p,pt);
      p[j] = pt;
      sum += pt;
    }
  }
}

/* ---------------------------------------------------------------------- */

void SolveGroup::update(int n, double *propensity)
{
  double pt = propensity[n];
  if (p[n] != pt) {
    if (p[n] > 0.0) num_active--;
    if (pt > 0.0) num_active++;
    sum -= p[n];
    groups->alter_element(n,p,pt);
    p[n] = pt;
    sum += pt;
  }
}

/* ---------------------------------------------------------------------- */

void SolveGroup::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}

/* ---------------------------------------------------------------------- */

int SolveGroup::event(double *pdt)
{
  int m;
  if (num_active == 0) {
    sum = 0.0;
    return -1;
  }

  m = groups->sample(p);
  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
