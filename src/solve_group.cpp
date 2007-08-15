/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_group.h"
#include "groups.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveGroup::SolveGroup(SPK *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg < 4) error->all("Illegal solve command");
  
  lo = atof(arg[1]);
  hi = atof(arg[2]);

  ngroups_in = 0;
  ngroups_flag = false;
  if (narg == 5) {
    ngroups_in = atoi(arg[3]); 
    ngroups_flag = true; 
    seed = atoi(arg[4]);
  }
  else seed = atoi(arg[3]);

  random = new RandomPark(seed);
  groups = new Groups(lo, hi, seed, ngroups_flag, ngroups_in);
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
  arg[3] = new char[16];
  arg[4] = new char[16];

  sprintf(arg[1],"%g",lo);
  sprintf(arg[2],"%g",hi);
  sprintf(arg[3],"%d",ngroups_in);
  sprintf(arg[4],"%d",seed);

  SolveGroup *ptr = new SolveGroup(spk,narg,arg);

  delete [] arg[1];
  delete [] arg[2];
  delete [] arg[3];
  delete [] arg[4];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveGroup::init(int n, double *propensity)
{
  nevents = n;
  nzeroes = 0;
  sum = 0.0;
  delete [] p;
  p = new double[n+10];

  last_size = n;

  for (int i = 0; i < n; i++) {
    double pt = propensity[i];
    if (lo > pt) {
      lo = pt; 
      fprintf(screen, "Lower bound violated. Reset to %g \n", lo);
    }
    else if (hi < pt)  {
      hi = pt; 
      fprintf(screen, "Upper bound violated. Reset to %g \n", hi);
    }
    p[i] = pt;
    sum += pt;
  }

  groups->partition_init(p,n,n+10);
}

/* ---------------------------------------------------------------------- */

void SolveGroup::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    int j = indices[i];
    double pt = p[j];
    if (propensity[j] != pt) {
      if (p[j] == 0.0) nzeroes--;
      if (propensity[j] == 0.0) nzeroes++;
      sum -= pt;
      groups->alter_element(j, p, pt);
      p[j] = propensity[j];
      sum += p[j];
    }
  }
}
/* ---------------------------------------------------------------------- */

void SolveGroup::update(int n, double *propensity)
{
  double pt = p[n];

  if (propensity[n] != pt) {
    if (pt == 0.0) nzeroes--;
    if (propensity[n] == 0.0) nzeroes++;
    sum -= pt;
    groups->alter_element(n, p, pt);
    p[n] = propensity[n];
    sum += p[n];
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

  if (nzeroes == nevents) return -1;

  m = groups->sample(p);
  *pdt = -1.0/sum * log(random->uniform());
  return m;
}

