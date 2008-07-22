/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
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

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPotts2d::AppPotts2d(SPPARKS *spk, int& narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  spinfile = NULL;

  // parse arguments

  if (narg < 5) error->all("Illegal app_style command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  int seed = atoi(arg[4]);
  random = new RandomPark(seed);
  init_style = RANDOM;

  // parse optional arguments
  // leave other arguments for child app

  int iarg = 5;
  int jarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"random") == 0) {
      init_style = RANDOM;
      iarg++;
    } else if (strcmp(arg[iarg],"read") == 0) {
      init_style = READ;
      iarg++;
      if (iarg >= narg) error->all("Illegal app_style command");
      int n = strlen(arg[iarg]) + 1;
      spinfile = new char[n];
      spinfile = strcpy(spinfile,arg[iarg]);
      iarg++;
    } else {
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
  delete random;
  delete [] spinfile;
}

