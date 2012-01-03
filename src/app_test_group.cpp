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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_test_group.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_mars.h"
#include "random_park.h"
#include "math.h"
#include "output.h"

using namespace SPPARKS_NS;

#define EPSILON 1.0e-10
//#define OUTPUT 1

/* ---------------------------------------------------------------------- */

AppTestGroup::AppTestGroup(SPPARKS *spk, int narg, char **arg) :
  App(spk, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal app_style command");

  nevents = atoi(arg[1]);
  ndep = atoi(arg[2]);
  pmax = atof(arg[3]);
  pmin = atof(arg[4]);
  tweak = atof(arg[5]);

  if (nevents == 0) error->all(FLERR,"Invalid event count for app_style test/group");
  if (pmin <= 0.0 || pmin >= pmax) 
    error->all(FLERR,"Invalid probability bounds for app_style test/group");
  if (tweak >= 100.0)
    error->all(FLERR,"Invalid probability delta for app_style test/group");
 
  pmax -= EPSILON*pmax;
  tweak = 2.0*tweak / 100.0;

  // optional args

  dep_graph = true;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lomem") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal app_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) dep_graph = false;
      else if (strcmp(arg[iarg+1],"no") == 0) dep_graph = true;
      else error->all(FLERR,"Illegal app_style command");
      iarg += 2;
    } else error->all(FLERR,"Illegal app_style command");
  }

  propensity = NULL;
  ndepends = NULL;
  depends = NULL;
  ran_dep = NULL;
  count = NULL;

  ncount = 0;

  // classes needed by this app

  random = new RandomPark(ranmaster->uniform());
}

/* ---------------------------------------------------------------------- */

AppTestGroup::~AppTestGroup()
{
  memory->destroy(propensity);
  memory->destroy(ndepends);
  memory->destroy(depends);
  memory->destroy(ran_dep);
  memory->destroy(count);
  delete random;
}

/* ---------------------------------------------------------------------- */

void AppTestGroup::input(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppTestGroup::init()
{
  // error check

  if (solve == NULL) error->all(FLERR,"No solver class defined");

  memory->destroy(ndepends);
  memory->destroy(depends);
  memory->destroy(ran_dep);

  if (dep_graph) build_dependency_graph();
  else memory->create(ran_dep,ndep,"test:ran_dep");

  memory->destroy(propensity);
  memory->create(propensity,nevents,"test:propensity");

  // allocate and zero event stats
  // initialize output

#ifdef OUTPUT
  memory->destroy(count);
  memory->create(count,nevents,"test:count");
  for (int m = 0; m < nevents; m++) count[m] = 0;

  output->init(time);
#endif
}

/* ---------------------------------------------------------------------- */

void AppTestGroup::setup()
{
  // compute initial propensity for each event

  double interval = log(pmax/pmin) / log(2.0);
  psum = 0.0;

  for (int m = 0; m < nevents; m++) {
    double p = pmax * pow(2.0,-random->uniform()*interval);
    p = MIN(p,pmax);
    p = MAX(p,pmin);
    propensity[m] = p;
    psum += propensity[m];
  }

  // initialize solver

  solve->init(nevents,propensity);

  // setup of output

#ifdef OUTPUT
  nextoutput = output->setup(time);
#endif
}

/* ----------------------------------------------------------------------
   iterate on solver
------------------------------------------------------------------------- */

void AppTestGroup::iterate()
{
  int i,m,n,ievent;
  double dt;

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    timer->stamp();
    ievent = solve->event(&dt);

#ifdef OUTPUT
    count[ievent]++;
#endif

    timer->stamp(TIME_SOLVE);

    propensity[ievent] = compute_propensity(ievent);
    solve->update(ievent,propensity);

    if (dep_graph) {
      for (m = 0; m < ndepends[ievent]; m++)
	propensity[depends[ievent][m]] = 
	  compute_propensity(depends[ievent][m]);
      solve->update(ndepends[ievent],depends[ievent],propensity);

    } else {
      n = static_cast<int>(ndep*random->uniform()) + 1;
      for (m = 0; m < n; m++) {
	i = static_cast<int> (nevents*random->uniform());
        ran_dep[m] = i;
	propensity[i] = compute_propensity(i);
      }
      solve->update(n,ran_dep,propensity);
    }

    timer->stamp(TIME_UPDATE);

    ncount++;
    time += dt;
    if (ncount >= stoptime) done = 1;
    else if (ievent < 0) done = 1;

#ifdef OUTPUT
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
#endif
  }

  timer->barrier_stop(TIME_LOOP);

  if (screen) fprintf(screen,"\nNumber of reactions, events = %d %d\n",
		      nevents,ncount);
  if (logfile) fprintf(logfile,"\nNumber of reactions, events = %d %d\n",
		       nevents,ncount);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppTestGroup::stats(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10g %10d",time,ncount);
  strpnt += strlen(strpnt);

  for (int m = 0; m < nevents; m++) {
    sprintf(strpnt," %d",count[m]);
    strpnt += strlen(strpnt);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppTestGroup::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s %s","Time","Step","Reaction-Counts");
}

/* ----------------------------------------------------------------------
   build random dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void AppTestGroup::build_dependency_graph()
{
  memory->create(ndepends,nevents,"test:ndepends");
  for (int m = 0; m < nevents; m++) 
    ndepends[m] = static_cast<int> (ndep*random->uniform()) + 1;

  int nmax = 1;
  for (int m = 0; m < nevents; m++)
    nmax = MAX(nmax,ndepends[m]);
  memory->create(depends,nevents,nmax,"test:depends");

  // set each dependency to random reaction other than self
  
  for (int m = 0; m < nevents; m++)
    for (int e = 0; e < ndepends[m]; e++) {
      depends[m][e] = m;
      while (depends[m][e] == m)
	depends[m][e] = static_cast<int> (nevents*random->uniform());
    }
}

/* ----------------------------------------------------------------------
   compute propensity of a single event
------------------------------------------------------------------------- */

double AppTestGroup::compute_propensity(int m)
{
  double p = propensity[m];
  p += p*tweak*(random->uniform()-0.5);
  
  p = MIN(p,pmax);
  p = MAX(p,pmin);
  return p;
}
