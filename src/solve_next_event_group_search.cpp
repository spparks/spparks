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

SolveNextEventGroupSearch::
SolveNextEventGroupSearch(SPK *spk, int narg, char **arg) : 
Solve(spk, narg, arg)
{
  if (narg < 4) error->all("Illegal solve command");
  
  lo = atof(arg[1]);
  hi = atof(arg[2]);

  ngroups_in = 0; ngroups_flag = false;
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

SolveNextEventGroupSearch::~SolveNextEventGroupSearch()
{
  delete random;
  delete groups;
  delete [] p;
}

/* ---------------------------------------------------------------------- */
//currently set to max size same as size
//to change, uncomment the commented lines
//and comment their complements
void SolveNextEventGroupSearch::init(int n, double *propensity)
{
  nevents = n;
  sum = 0;
  delete [] p;
  //p = new double[2*n];
  p = new double[n+10];
  //p = propensity;

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

  //  groups->partition_init(propensity, n, 2*n);
  groups->partition_init(p, n, n+10);
  //groups->partition_init(p, n, n+10);
}

/* ---------------------------------------------------------------------- */

void SolveNextEventGroupSearch::update(int n, int *indices, double *propensity)
//void SolveNextEventGroupSearch::update(int n, int *indices, double *pold)
{
  for (int i = 0; i < n; i++) {
    int j = indices[i];
    //    update(j, propensity);
    // double pt = pold[i];
    double pt = p[j];
    if(propensity[j]!=pt){
      sum -= pt;
      groups->alter_element(j, p, pt);
      p[j] = propensity[j];
      sum +=  p[j];
    }
  }
}
/* ---------------------------------------------------------------------- */

void SolveNextEventGroupSearch::update(int n, double *propensity)
//void SolveNextEventGroupSearch::update(int n, double *pold)
{
  double pt = p[n];
  //double pt = pold[0];

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
  if(propensity[n]!=pt){
    //sum -= p[n];
    sum -= pt;
    //groups->alter_element(n, p, pt);
    groups->alter_element(n, p, pt);
    p[n] = propensity[n];
    //sum +=  pt;
    sum += p[n];
  }
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

