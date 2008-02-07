/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_group2.h"
#include "groups2.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SolveGroup2::SolveGroup2(SPK *spk, int narg, char **arg) :
  Solve(spk, narg, arg)
{
  if (narg < 4) error->all("Illegal solve command");
  
  lo = atof(arg[1]);
  hi = atof(arg[2]);

  ngroups_in = 0;
  ngroups_flag = false;
  if (narg == 5) {
    ngroups_in = atoi(arg[3]);
    if (ngroups_in) ngroups_flag = true;
    seed = atoi(arg[4]);
  }
  else seed = atoi(arg[3]);

  random = new RandomPark(seed);
  groups = new Groups2(spk,lo,hi,seed,ngroups_flag,ngroups_in);
  p = NULL;
}

/* ---------------------------------------------------------------------- */

SolveGroup2::~SolveGroup2()
{
  delete random;
  delete groups;
  delete [] p;
}

/* ---------------------------------------------------------------------- */

SolveGroup2 *SolveGroup2::clone()
{
  int narg = 5;
  char *arg[5];

  arg[0] = style;
  arg[1] = new char[16];
  arg[2] = new char[16];
  arg[3] = new char[16];
  arg[4] = new char[16];

  sprintf(arg[1],"%g",lo);
  sprintf(arg[2],"%g",hi);
  sprintf(arg[3],"%d",ngroups_in);
  sprintf(arg[4],"%d",seed);

  SolveGroup2 *ptr = new SolveGroup2(spk,narg,arg);

  delete [] arg[1];
  delete [] arg[2];
  delete [] arg[3];
  delete [] arg[4];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveGroup2::init(int n, double *propensity)
{
  nevents = n;
  num_active = 0;

  delete [] p;
  p = new double[n];

  sum = 0.0;
  for (int i = 0; i < n; i++) {
    p[i] = MAX(propensity[i],lo);
    p[i] = MIN(p[i],hi);
    sum += p[i];
  }

  groups->partition_init(p,n);
}

/* ---------------------------------------------------------------------- */

void SolveGroup2::update(int n, int *indices, double *propensity)
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

void SolveGroup2::update(int n, double *propensity)
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

void SolveGroup2::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}

/* ---------------------------------------------------------------------- */

int SolveGroup2::event(double *pdt)
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
