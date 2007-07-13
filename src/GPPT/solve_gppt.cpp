
/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_gppt.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"
#include <iostream>

using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveGPPT::SolveGPPT(SPK *spk, int narg, char **arg) : Solve(spk, narg, arg) 
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);

  //  cout<< "Using solver GPPT."<<endl;
}

/* ---------------------------------------------------------------------- */

SolveGPPT::~SolveGPPT()
{
   delete random;
}

/* ---------------------------------------------------------------------- */

void SolveGPPT::resize(int new_size, double *propensity)
{
  init(new_size, propensity);
}
/* ---------------------------------------------------------------------- */

int SolveGPPT::event(double *pdt)
{
  int m = 0;
  return m;
}
