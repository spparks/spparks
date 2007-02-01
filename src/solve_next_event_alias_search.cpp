/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_next_event_alias_search.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveNextEventAliasSearch::SolveNextEventAliasSearch
(SPK *spk, int narg, char **arg) : Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
}

/* ---------------------------------------------------------------------- */

SolveNextEventAliasSearch::~SolveNextEventAliasSearch()
{
  delete random;
  free_arrays();
}

/* ---------------------------------------------------------------------- */

void SolveNextEventAliasSearch::init(int n, double *propensity)
{
  int i;

  if(allocated) free_arrays();
  allocated = 1;
  nevents = n;

  j = new int[n];
  p = new double[n];
  q = new double[n];
  hilo = new int[n];

  build_alias_table(nevents, propensity);
}
/* ----------------------------------------------------------------------
   build table
------------------------------------------------------------------------- */

void SolveNextEventAliasSearch::build_alias_table(int n, double *propensity)
{
  int i, m, k;

  sum = 0.0;

  for (i = 0; i < n; i++) sum += propensity[i];

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
int SolveNextEventAliasSearch::sort_hilo()
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
void SolveNextEventAliasSearch::check_table_consistency()
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
void SolveNextEventAliasSearch::table_dump(int n)
{
  int i;
  fprintf(screen,"----- table dump: -----\n");
  
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
  fprintf(screen,"-----------------------\n");
}

/* ---------------------------------------------------------------------- */

void SolveNextEventAliasSearch::update(int n, int *indices, double *propensity)
{

  //remember to update sum

  build_alias_table(nevents, propensity);
}

/* ---------------------------------------------------------------------- */

int SolveNextEventAliasSearch::event(double *pdt)
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
void SolveNextEventAliasSearch::free_arrays()
{

  delete [] p;
  delete [] q;
  delete [] j;
}
