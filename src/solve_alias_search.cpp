/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_alias_search.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

#include <iostream>

using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveAliasSearch::SolveAliasSearch(SPK *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
  p = NULL;
  q = NULL;
  j = NULL;
}

/* ---------------------------------------------------------------------- */

SolveAliasSearch::~SolveAliasSearch()
{
  delete random;
  free_arrays();
}

/* ---------------------------------------------------------------------- */

void SolveAliasSearch::init(int n, double *propensity)
{
  int i;

  if(allocated) free_arrays();
  allocated = 1;
  nevents = n;

  prob = new double[n];
  j = new int[n];
  p = new double[n];
  q = new double[n];
  hilo = new int[n];

  sum = 0.0;
  for (i = 0; i < n; i++) {sum += propensity[i]; prob[i] = propensity[i];}

  build_alias_table(nevents, propensity);
  //  table_dump(nevents);
}

/* ----------------------------------------------------------------------
   build table
------------------------------------------------------------------------- */

void SolveAliasSearch::build_alias_table(int n, double *propensity)
{
  int i, m, k;


  for (i = 0; i < n; i++){
    p[i] = propensity[i]/sum;
    q[i] = (double)n*p[i];
    j[i] = i;
    hilo[i]=i;
  }
  sk = sort_hilo();
  sl = sk - 1;

  while (sl >  -1 & sk < nevents){

    k = hilo[sk];
    l = hilo[sl];

    j[l] = k;
    q[k] += q[l] - 1.0;

    if(q[k]<1.0){
      hilo[sl] = hilo[sk];
      sk++;
    }
    else  sl--;

    //    table_dump(nevents);   //diagnostic
  }
  //  check_table_consistency();  //diagnostic
}

/* ---------------------------------------------------------------------- */
//sort hilo array
/* ---------------------------------------------------------------------- */

int SolveAliasSearch::sort_hilo()
{
  int i, m;
  int tmp;

  for (i=0;i<nevents-1;i++)
    for (m=i+1;m<nevents;m++)
      if(q[hilo[i]]>q[hilo[m]]) {
	tmp = hilo[i];
	hilo[i] = hilo[m];
	hilo[m] = tmp;
      }
  i=0;
  while (i< nevents && q[hilo[i]]<1.0) i++;

  return i;
}

/* ---------------------------------------------------------------------- */
//diagnostic table check
/* ---------------------------------------------------------------------- */

void SolveAliasSearch::check_table_consistency()
{
  double psum;

  for (int i=0;i<nevents;i++){
    psum = 0.0;
    psum += q[i];
    for (k=0;k<nevents;k++)
      if(j[k]==i) psum += 1.0 - q[k];
    if ((double)nevents*p[i]-psum >1.0e-12) 
      fprintf(screen,"propensity %d incomplete in table \n",i);
  }
}

/* ---------------------------------------------------------------------- */
//diagnostic table dump
/* ---------------------------------------------------------------------- */

void SolveAliasSearch::table_dump(int n)
{
  int i;
  fprintf(screen,
	  "%%%%%%%%%%%%%%%%%%%%%%%%%%%% table dump: %%%%%%%%%%%%%%%%%%%\n");
  
  fprintf(screen,"i: ");
  for (i = 0; i < nevents; i++)fprintf(screen," %d ",i);
  fprintf(screen,"\n");
  fprintf(screen,"j: ");
  for (i = 0; i < nevents; i++)fprintf(screen," %d ",j[i]);
  fprintf(screen,"\n");
  fprintf(screen,"p: ");
  for (i = 0; i < nevents; i++)fprintf(screen," %f ",p[i]);
  fprintf(screen,"\n");
  fprintf(screen,"q: ");
  for (i = 0; i < nevents; i++)fprintf(screen," %f ",q[i]);
  fprintf(screen,"\n");
  fprintf(screen,
	  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
}

/* ---------------------------------------------------------------------- */

void SolveAliasSearch::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    int current_index = indices[i];
    sum -= prob[current_index];
    prob[current_index] = propensity[current_index];
    sum +=  propensity[current_index];
  }

  if (n > 0) {build_alias_table(nevents, propensity);}
}

/* ---------------------------------------------------------------------- */

void SolveAliasSearch::update(int n, double *propensity)
{
  sum -= prob[n];
  prob[n] = propensity[n];
  sum +=  propensity[n];
  
  build_alias_table(nevents, propensity);
}

/* ---------------------------------------------------------------------- */

void SolveAliasSearch::resize(int new_size, double *propensity)
{
  init(new_size, propensity);
}

/* ---------------------------------------------------------------------- */

int SolveAliasSearch::event(double *pdt)
{
  int i;

  *pdt = -1.0/sum * log(random->uniform());

  i = (int)((double)nevents * random->uniform());

  //  fprintf(screen,"event \n");

  if(random->uniform()<q[i]) return i;
  
  return j[i];
}

/* ----------------------------------------------------------------------
   free arrays used by solver
------------------------------------------------------------------------- */
void SolveAliasSearch::free_arrays()
{

  delete [] p;
  delete [] q;
  delete [] j;
}
