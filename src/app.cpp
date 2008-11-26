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

#include "string.h"
#include "app.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

App::App(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
  appclass = GENERAL;
  time_eps = 1.0e-6;
}

/* ---------------------------------------------------------------------- */

App::~App()
{
  delete [] style;
}
