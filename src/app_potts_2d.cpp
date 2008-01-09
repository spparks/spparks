/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_2d.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppPotts2d::AppPotts2d(SPK *spk, int& narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{

  spinfile = NULL;

  // parse arguments

  if (narg < 5) error->all("Invalid app_style potts/2d command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  seed = atoi(arg[4]);
  random = new RandomPark(seed);
  init_style = RANDOM;

  // parse optional arguments

  int iarg = 5;
  int jarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"random") == 0) {
      init_style = RANDOM;
      iarg ++;
    } else if (strcmp(arg[iarg],"read") == 0) {
      init_style = READ;
      iarg ++;
      if (iarg >= narg) error->all("Illegal app_style potts/2d command");
      int n = strlen(arg[iarg]) + 1;
      spinfile = new char[n];
      spinfile = strcpy(spinfile,arg[iarg]);
      iarg ++;
    } else {
      // leave other arguments for child app
      arg[jarg] = arg[iarg];
      iarg++;
      jarg++;
    }
  }
  narg = jarg;
}

/* ---------------------------------------------------------------------- */

AppPotts2d::~AppPotts2d()
{
  delete [] spinfile;
}

