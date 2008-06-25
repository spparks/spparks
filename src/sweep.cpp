/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "string.h"
#include "sweep.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Sweep::Sweep(SPPARKS *spk, int narg, char **arg) : SysPtr(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
}

/* ---------------------------------------------------------------------- */

Sweep::~Sweep()
{
  delete [] style;
}
