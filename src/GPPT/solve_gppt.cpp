
/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_gppt.h"
#include "random_park.h"
#include "error.h"

using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveGPPT::SolveGPPT(SPK *spk, int narg, char **arg) : Solve(spk, narg, arg) 
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
}

/* ---------------------------------------------------------------------- */

SolveGPPT::~SolveGPPT()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

SolveGPPT *SolveGPPT::clone()
{
  int narg = 2;
  char *arg[2];
  arg[0] = style;
  arg[1] = new char[16];
  sprintf(arg[1],"%d",seed);

  SolveGPPT *ptr = new SolveGPPT(spk,narg,arg);

  delete [] arg[1];
  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveGPPT::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}
/* ---------------------------------------------------------------------- */

int SolveGPPT::event(double *pdt)
{
  int m = 0;
  return m;
}
