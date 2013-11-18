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
#include "domain.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

SolveGroup::SolveGroup(SPPARKS *spk, int narg, char **arg) :
  Solve(spk, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal solve command");
  
  hi = atof(arg[1]);
  lo = atof(arg[2]);

  if (lo <= 0.0 || lo >= hi)
    error->all(FLERR,"Invalid probability bounds for solve_style group");

  ngroups = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ngroup") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal solve_style group command");
      ngroups = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal solve_style group command");
  }

  // each proc uses different initial RNG seed

  random = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  random->reset(seed,spk->domain->me,100);

  groups = new Groups(spk,hi,lo,ngroups);
  p = NULL;

  nroundlo = nroundhi = 0;
  lomax = lo;
  himax = hi;
}

/* ---------------------------------------------------------------------- */

SolveGroup::~SolveGroup()
{
  round_check();

  delete random;
  delete groups;
  memory->destroy(p);
}

/* ---------------------------------------------------------------------- */

SolveGroup *SolveGroup::clone()
{
  int narg = 5;
  char *arg[5];

  arg[0] = style;
  arg[1] = new char[16];
  arg[2] = new char[16];
  arg[3] = (char *) "ngroup";
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

/* ----------------------------------------------------------------------
   bound all non-zero propensities between LO and HI inclusive
   zero propensities are not changed and do not contribute to num_active
------------------------------------------------------------------------- */

void SolveGroup::init(int n, double *propensity)
{
  nevents = n;
  num_active = 0;

  memory->destroy(p);
  memory->create(p,n,"solve/group:p");

  sum = 0.0;
  for (int i = 0; i < n; i++) {
    if (propensity[i] > 0.0) {
      num_active++;
      p[i] = propensity[i];
      if (p[i] < lo) {
	nroundlo++;
	lomax = MIN(lomax,p[i]);
	p[i] = lo;
      }
      if (p[i] > hi) {
	nroundhi++;
	himax = MAX(himax,p[i]);
	p[i] = hi;
      }
      sum += p[i];
    } else p[i] = 0.0;
  }

  groups->partition(p,n);

  round_check();
  nroundlo = nroundhi = 0;
  lomax = lo;
  himax = hi;
}

/* ----------------------------------------------------------------------
   update a list of propensity values
   invoke alter_element() only if a propensity changes
   bound all non-zero propensities between LO and HI inclusive
   adjust num_active based on old p[j] vs new pt propensity values
------------------------------------------------------------------------- */

void SolveGroup::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    int j = indices[i];
    double pt = propensity[j];
    if (p[j] != pt) {
      if (p[j] == 0.0) num_active++;
      if (pt == 0.0) num_active--;
      else {
	if (pt < lo) {
	  nroundlo++;
	  lomax = MIN(lomax,pt);
	  pt = lo;
	}
	if (pt > hi) {
	  nroundhi++;
	  himax = MIN(himax,pt);
	  pt = hi;
	}
      }
      sum -= p[j];
      groups->alter_element(j,p,pt);
      p[j] = pt;
      sum += pt;
    }
  }
}

/* ----------------------------------------------------------------------
   update a single propensity value
   invoke alter_element() only if propensity changes
   bound non-zero propensity between LO and HI inclusive
   adjust num_active based on old p[j] vs new pt propensity value
------------------------------------------------------------------------- */

void SolveGroup::update(int n, double *propensity)
{
  double pt = propensity[n];
  if (p[n] != pt) {
    if (p[n] == 0.0) num_active++;
    if (pt == 0.0) num_active--;
    else {
      if (pt < lo) {
	nroundlo++;
	lomax = MIN(lomax,pt);
	pt = lo;
      }
      if (pt > hi) {
	nroundhi++;
	himax = MIN(himax,pt);
	pt = hi;
      }
    }
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

/* ---------------------------------------------------------------------- */

void SolveGroup::round_check()
{
  int nlo,nhi;
  double lomaxall,himaxall;
  MPI_Allreduce(&nroundlo,&nlo,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nroundhi,&nhi,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&lomax,&lomaxall,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&himax,&himaxall,1,MPI_DOUBLE,MPI_MAX,world);

  int me;
  MPI_Comm_rank(world,&me);
  if ((nlo || nhi) && me == 0) {
    char str[128];
    sprintf(str,"%d propensities were reset to lo value, max lo = %g",
	    nlo,lomaxall);
    error->warning(FLERR,str);
    sprintf(str,"%d propensities were reset to hi value, max hi = %g",
	    nhi,himaxall);
    error->warning(FLERR,str);
  }
}
