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
}

/* ---------------------------------------------------------------------- */

Solve::~Solve()
{
  delete [] style;
}
