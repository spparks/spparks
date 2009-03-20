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
#include "stdlib.h"
#include "app.h"
#include "timer.h"
#include "finish.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

App::App(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  appclass = GENERAL;
  time = 0.0;
  first_run = 1;
}

/* ---------------------------------------------------------------------- */

App::~App()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

void App::run(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal run command");

  stoptime = time + atof(arg[0]);
  if (stoptime < time) error->all("Illegal run command");

  // read optional args

  int uptoflag = 0;
  int preflag = 1;
  int postflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) error->all("Illegal run command");
      uptoflag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) preflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) preflag = 1;
      else error->all("Illegal run command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) postflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) postflag = 1;
      else error->all("Illegal run command");
      iarg += 2;
    } else error->all("Illegal run command");
  }

  // adjust stoptime if upto was specified

  if (uptoflag) stoptime -= time;
  if (uptoflag && stoptime < 0.0)
    error->all("Run upto value is before current time");

  // perform a single run via app's init(), setup(), and iterate()
  // if pre or 1st run, do app init
  // setup computes initial propensities
  // if post, do full Finish, else just print time

  if (preflag || first_run) {
    init();
    first_run = 0;
  }
  timer->init();
  setup();

  if (stoptime > time) iterate();

  Finish finish(spk,postflag);
}
