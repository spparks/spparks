/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_alias.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveAlias::SolveAlias(SPK *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  p = NULL;
  q = NULL;
  j = NULL;
}

/* ---------------------------------------------------------------------- */

SolveAlias::~SolveAlias()
{
  delete random;
  free_arrays();
}

/* ---------------------------------------------------------------------- */

SolveAlias *SolveAlias::clone()
{
  int narg = 2;
  char *arg[2];
  arg[0] = style;
  arg[1] = new char[16];
  sprintf(arg[1],"%d",seed);

  SolveAlias *ptr = new SolveAlias(spk,narg,arg);

  delete [] arg[1];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveAlias::init(int n, double *propensity)
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

  num_active = 0;
  sum = 0.0;
  for (i = 0; i < n; i++) {
    if (propensity[i] > 0.0) num_active++;
    prob[i] = propensity[i];
    sum += propensity[i];
  }

  build_alias_table(nevents, propensity);
  //  table_dump(nevents);
}

/* ----------------------------------------------------------------------
   build table
------------------------------------------------------------------------- */

void SolveAlias::build_alias_table(int n, double *propensity)
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

int SolveAlias::sort_hilo()
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

void SolveAlias::check_table_consistency()
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

void SolveAlias::table_dump(int n)
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

void SolveAlias::update(int n, int *indices, double *propensity)
{
  int m;
  for (int i = 0; i < n; i++) {
    m = indices[i];
    if (prob[m] > 0.0) num_active--;
    if (propensity[m] > 0.0) num_active++;
    sum -= prob[m];
    prob[m] = propensity[m];
    sum +=  propensity[m];
  }
  if (n > 0) build_alias_table(nevents, propensity);
}

/* ---------------------------------------------------------------------- */

void SolveAlias::update(int n, double *propensity)
{
  if (prob[n] > 0.0) num_active--;
  if (propensity[n] > 0.0) num_active++;
  sum -= prob[n];
  prob[n] = propensity[n];
  sum += propensity[n];
  
  build_alias_table(nevents, propensity);
}

/* ---------------------------------------------------------------------- */

void SolveAlias::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}

/* ---------------------------------------------------------------------- */

int SolveAlias::event(double *pdt)
{
  int i;

  if (num_active == 0) return -1;

  *pdt = -1.0/sum * log(random->uniform());

  i = (int)((double)nevents * random->uniform());
  if (random->uniform() < q[i]) return i;
  return j[i];
}

/* ----------------------------------------------------------------------
   free arrays used by solver
------------------------------------------------------------------------- */
void SolveAlias::free_arrays()
{

  delete [] p;
  delete [] q;
  delete [] j;
}

