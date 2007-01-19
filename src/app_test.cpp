/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_test.h"
#include "spk.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"

using namespace SPPARKS;


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppTest::AppTest(SPK *spk, int narg, char **arg) : App(spk, narg, arg)
{
  if (narg != 1) error->all("Invalid app_style test command");

  ndepends = NULL;
  depends = NULL;
  propensity = NULL;

  nevents = 0;
  n_event_types = 0;
  ntimestep = 0;
  time = 0.0;
  stoptime = 0.0;
  // classes needed by this app
  int seed = 1234567;
  random = new RandomPark(seed);

  delete timer;
  timer = new Timer(spk);
}

/* ---------------------------------------------------------------------- */

AppTest::~AppTest()
{
  delete random;
  delete [] propensity;
  delete [] ndepends;
  memory->destroy_2d_int_array(depends);
}

/* ---------------------------------------------------------------------- */

void AppTest::init()
{

  if (n_event_types == 0)
    error->all("No events defined for test app");
  if (nevents == 0)
    error->all("Zero events defined for test app");
  // determine event dependencies

  delete [] ndepends;
  memory->destroy_2d_int_array(depends);
  ndepends = new int[nevents];
  build_dependency_graph();

  // compute initial propensity for each event
  // inform Nfold solver

  delete [] propensity;
  propensity = new double[nevents];
  for (int m = 0; m < nevents; m++) propensity[m] = compute_propensity(m);
  solve->init(nevents,propensity);

  // allocate and zero event stats
  count = new int[nevents];
  for (int t = 0; t< nevents; t++)count[t] = 0;


  // print stats header

  if (screen) {
    fprintf(screen,"Step Time Counts");
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"Step Time Counts");
    fprintf(logfile,"\n");
  }
  stats();

  // setup future calls to stats()

  stats_time = time + stats_delta;
  if (stats_delta == 0.0) stats_time = stoptime;
}

/* ---------------------------------------------------------------------- */

void AppTest::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  else if (strcmp(command,"event") == 0) set_event(narg,arg);
  else if (strcmp(command,"run") == 0) run(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else error->all("Invalid command");
}

/* ----------------------------------------------------------------------
   perform a run
------------------------------------------------------------------------- */

void AppTest::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

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

void AppTest::iterate()
{
  int d,ievent;
  double dt;

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    ntimestep++;

    timer->stamp();
    ievent = solve->event(&dt);
    count[ievent] ++;
    timer->stamp(TIME_SOLVE);

    
    // update propensity table
    // inform solver of changes

    // first test simply replacing an old propensity of an event by a new
    // this is equivalent to removing one event and adding another
    // so far the dependency is only self

    for (d = 0; d < ndepends[ievent]; d++)
      propensity[depends[ievent][d]] = compute_propensity(depends[ievent][d]);

    solve->update(ndepends[ievent],depends[ievent],propensity);

    // update time by dt

    time += dt;
    if (time >= stoptime) done = 1;
    else if (ievent < 0) done = 1;

    if (time > stats_time || done) {
      timer->stamp();
      stats();
      stats_time += stats_delta;
      timer->stamp(TIME_OUTPUT);
    }
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppTest::stats()
{
  int ssum = 0;
  int i;

  for (i = 0; i< nevents; i++) ssum += count[i];

  if (screen) {
    fprintf(screen,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(screen,"%g ", 
	      (double)count[i]/(double)ssum);
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(logfile,"%g ", 
	      (double)count[i]/(double)ssum);
    fprintf(logfile,"\n");
    
  }
}
/* ---------------------------------------------------------------------- */

void AppTest::set_event(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal event command");

  for(int iarg = 0; iarg < narg; iarg++)
    fprintf(screen,"number of events: %s  ",arg[iarg]);
    fprintf(screen,"\n");

    n_event_types ++;
    nevents = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppTest::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
}

/* ----------------------------------------------------------------------
   build dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void AppTest::build_dependency_graph()
{
  int m;

   int nmax = 1;
   for (m = 0; m < nevents; m++)
     nmax = MAX(nmax,ndepends[m]);

   depends = memory->create_2d_int_array(nevents,nmax,
					      "test:depends");
   // zero the dependencies
   // include self
   for (m = 0; m < nevents; m++) {
     ndepends[m] = 1;
     depends[m][0]= m;
   }

}

/* ----------------------------------------------------------------------
   compute propensity of a single event
------------------------------------------------------------------------- */

double AppTest::compute_propensity(int m)
{
  //uniform
  //  double p=.1;
  //random uniform
  //double p=random->uniform();
  //linear
  double p = (double)(m+1);

  return p;
}
