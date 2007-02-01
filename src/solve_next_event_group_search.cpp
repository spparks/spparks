/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_next_event_group_search.h"
#include "groups.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveNextEventGroupSearch::SolveNextEventGroupSearch(SPK *spk, int narg, char **arg) : Solve(spk, narg, arg)
{
  if (narg != 4) error->all("Illegal solve command");
  
  lo = atof(arg[1]);
  hi = atof(arg[2]);
  seed = atoi(arg[3]);

  random = new RandomPark(seed);
  groups = new Groups(lo, hi, seed);

  fprintf(screen,"Using groups search algorithm to generate events.\n");
}

/* ---------------------------------------------------------------------- */

SolveNextEventGroupSearch::~SolveNextEventGroupSearch()
{
  delete random;
  delete groups;
  delete [] p;
}

/* ---------------------------------------------------------------------- */

void SolveNextEventGroupSearch::init(int n, double *propensity)
{
  nevents = n;
  sum = 0;
  if (p!=NULL) delete [] p;
  p = new double[2*n];
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

  groups->partition_init(propensity, n, 2*n);
}

/* ---------------------------------------------------------------------- */

void SolveNextEventGroupSearch::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    int j = indices[i];
    //    update(j, propensity);
    double pt = propensity[j];
    sum -= p[j];
    groups->alter_element(j, p, pt);
    p[j] = pt;
    sum +=  pt;
  }
}
/* ---------------------------------------------------------------------- */

void SolveNextEventGroupSearch::update(int n, double *propensity)
{
    double pt = propensity[n];
//     if (lo > pt) {
//       lo = pt; 
//       fprintf(screen, "Lower bound violated. Reset to %g \n", lo);
//       init(n, propensity);
//     }
//     else if (hi < pt)  {
//       hi = pt; 
//       fprintf(screen, "Upper bound violated. Reset to %g \n", hi);
//       init(n, propensity);
//     }
//     else{
      sum -= p[n];
      groups->alter_element(n, p, pt);
      p[n] = pt;
      sum +=  pt;
      //    }
}
/* ---------------------------------------------------------------------- */


void SolveNextEventGroupSearch::resize(int new_size, double *propensity)
{
  init(new_size, propensity);
}
/* ---------------------------------------------------------------------- */
int SolveNextEventGroupSearch::event(double *pdt)
{
  int m;

  if (sum == 0.0) return -1;
  m = groups->sample(p);

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
