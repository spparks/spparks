/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_test2_group.h"
#include "spk.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "math.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppTest2Group::AppTest2Group(SPK *spk, int narg, char **arg) :
  App(spk, narg, arg)
{
  if (narg != 1) error->all("Invalid app_style test2/group command");

  propensity = NULL;
  ndepends = NULL;
  depends = NULL;
  ran_dep = NULL;

  nevents = 0;
  dep_graph_flag = true;
  time = 0.0;

  // classes needed by this app

  int seed = 123124;
  random = new RandomPark(seed);
}

/* ---------------------------------------------------------------------- */

AppTest2Group::~AppTest2Group()
{
  delete [] propensity;
  delete [] ndepends;
  memory->destroy_2d_int_array(depends);
  delete [] ran_dep;

  delete random;
}

/* ---------------------------------------------------------------------- */

void AppTest2Group::init()
{
  if (nevents == 0)
    error->all("Zero events defined for test app");

  delete [] ndepends;
  memory->destroy_2d_int_array(depends);
  delete [] ran_dep;

  if (dep_graph_flag) build_dependency_graph();
  else ran_dep = new int[ndep];

  // compute initial propensity for each event
  // inform solver

  delete [] propensity;
  propensity = new double[nevents];

  psum = 0;
  for (int m = 0; m < nevents; m++) {
    double tp = pow(2.0,-random->uniform()*20);
    propensity[m] = MAX(tp,1.0e-6);
    psum += propensity[m];
  }

  // print header

  if (screen) {
    fprintf(screen,"Starting Test");
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"Starting Test");
    fprintf(logfile,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void AppTest2Group::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  else if (strcmp(command,"event") == 0) set_event(narg,arg);
  else if (strcmp(command,"run") == 0) run(narg,arg);
  else error->all("Invalid command");
}

/* ----------------------------------------------------------------------
   perform a run
------------------------------------------------------------------------- */

void AppTest2Group::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  nlimit = atoi(arg[0]);

  // error check

  if (solve == NULL) error->all("No solver class defined");

  // init classes used by this app
  
  init();
  solve->init(nevents,propensity);
  timer->init();

  // perform the run

  iterate();

  // final statistics

  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on solver
------------------------------------------------------------------------- */

void AppTest2Group::iterate()
{
  int i,m,n,ievent,rdp,tdep;
  double dt;

  // comment out if don't want event counts

  //count = new int[nevents];
  //for (i = 0; i < nevents; i++) count[i] = 0;

  ntimestep = 0;
  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    ntimestep++;

    timer->stamp();
    ievent = solve->event(&dt);

    // comment out if don't want event counts
    //count[ievent]++;

    timer->stamp(TIME_SOLVE);

    propensity[ievent] = compute_propensity(ievent);
    solve->update(ievent,propensity);

    if (dep_graph_flag) {
      for (m = 0; m < ndepends[ievent]; m++)
	propensity[depends[ievent][m]] = 
	  compute_propensity(depends[ievent][m]);
      solve->update(ndepends[ievent],depends[ievent],propensity);

    } else {
      n = static_cast<int>(ndep*random->uniform()) + 1;
      for (m = 0; m < n; m++) {
	i = static_cast<int>(nevents*random->uniform());
        ran_dep[m] = i;
	propensity[i] = compute_propensity(i);
      }
      solve->update(n,ran_dep,propensity);
    }

    timer->stamp(TIME_UPDATE);

    time += dt;
    if (ntimestep >= nlimit) done = 1;
    else if (ievent < 0) done = 1;
  }

  timer->barrier_stop(TIME_LOOP);

  // comment out if don't want event counts

  //stats();
  //delete [] count;

  if (screen) fprintf(screen,"\nNumber of reactions, events = %d %d\n",
		      nevents,ntimestep);
  if (logfile) fprintf(logfile,"\nNumber of reactions, events = %d %d\n",
		       nevents,ntimestep);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppTest2Group::stats()
{
  if (screen) {
    fprintf(screen,"Reactions counts:\n");
    for (int m = 0; m < nevents; m++) fprintf(screen," %d",count[m]);
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"Reactions counts:\n");
    for (int m = 0; m < nevents; m++) fprintf(logfile," %d",count[m]);
    fprintf(logfile,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void AppTest2Group::set_event(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal event command");

  if (narg > 4)
    if (strcmp(arg[3],"lo_mem")==0) {
      dep_graph_flag = false;
      random->init(atoi(arg[4]));
    } else error->all("Illegal event command");

  nevents = atoi(arg[0]);
  ndep = atoi(arg[1]);
  tweak = atof(arg[2]);
  if (tweak > 0) tweak /= 100.0;
}

/* ----------------------------------------------------------------------
   build random dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void AppTest2Group::build_dependency_graph()
{
  ndepends = new int[nevents];
  for (int m = 0; m < nevents; m++) 
    ndepends[m] = static_cast<int>(ndep*random->uniform()) + 1;

  int nmax = 1;
  for (int m = 0; m < nevents; m++)
    nmax = MAX(nmax,ndepends[m]);
  depends = memory->create_2d_int_array(nevents,nmax,"test:depends");

  // set each dependency to random reaction other than self
  
  for (int m = 0; m < nevents; m++)
    for (int e = 0; e < ndepends[m]; e++) {
      depends[m][e] = m;
      while (depends[m][e] == m)
	depends[m][e] = static_cast<int>(nevents*random->uniform());
    }
}

/* ----------------------------------------------------------------------
   compute propensity of a single event
------------------------------------------------------------------------- */

double AppTest2Group::compute_propensity(int m)
{
  double p = propensity[m];
  p += p*tweak*(random->uniform()-0.5);
  
  if (p > 1.0) p = 0.99999;
  else if (p < 1.0e-6) p = 1.0e-6;
  return p;
}
