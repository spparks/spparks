/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "string.h"
#include "solve.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

Solve::Solve(SPK *spk, int narg, char **arg) : SysPtr(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
  sum = 0.0;
}

/* ---------------------------------------------------------------------- */

Solve::~Solve()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

double Solve::get_total_propensity()
{
  return sum;
}

/* ---------------------------------------------------------------------- */

int Solve::get_num_active()
{
  return num_active;
}
