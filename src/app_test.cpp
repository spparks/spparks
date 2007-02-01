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
  count = NULL;

  nevents = 0;
  n_event_types = 0;
  ntimestep = 0;
  time = 0.0;
  stoptime = 0.0;
  // classes needed by this app
  int seed = 123124;
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
  delete [] count;
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
  psum = 0;
  delete [] propensity;
  propensity = new double[nevents];
  for (int m = 0; m < nevents; m++) {
    propensity[m] = compute_propensity(m);
    psum += propensity[m];
  }
  //  solve->init(nevents,propensity);
  
  // allocate and zero event stats
  if(count != NULL) delete [] count;
  count = new int[2*nevents];
  for (int t = 0; t< 2*nevents; t++) count[t] = 0;

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
  // uncomment to control number of total events
  int nev = 0;
  double dt;

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    //uncomment to control total number of events
    nev ++;

    ntimestep++;
    timer->stamp();
    ievent = solve->event(&dt);
    count[ievent] ++;
    timer->stamp(TIME_SOLVE);

    // update propensity table
    // inform solver of changes

    //update event propensity
    propensity[ievent] = compute_propensity(ievent);
    solve->update(ievent, propensity);
    
    //update dependencies
    for (d = 0; d < ndepends[ievent]; d++)
      propensity[depends[ievent][d]] = compute_propensity(depends[ievent][d]);
    
    solve->update(ndepends[ievent],depends[ievent],propensity);
    
    timer->stamp(TIME_UPDATE);
    // update time by dt

    time += dt;
    if (time >= stoptime) done = 1;
    else if (ievent < 0) done = 1;
    // uncomment to control total number of events
    else if(nev > 100000) done = 1;

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
  int i;
  ssum = 0;
  double deviation;
  double max_deviation = 0.0;
  int max_i;

  for (i = 0; i< nevents; i++) ssum += count[i];

  if (screen)
    //    fprintf(screen,"Deviations: \n");
  if (screen) {
    //    fprintf(screen,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++){
      deviation = ((double)count[i]/(double)ssum-propensity[i]/psum)/
	((double)count[i]/(double)ssum+propensity[i]/psum);

      if(deviation > max_deviation){
	max_deviation = deviation;
	max_i = i;
      }
      //      fprintf(screen,"%g ", deviation);
    }
    //    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(logfile,"%g ", 
	      (double)count[i]/(double)ssum-propensity[i]/psum);
    fprintf(logfile,"\n");
    
  }
  if(max_i>nevents-1) max_i = nevents - 1;
  if (screen)
    fprintf(screen,"Maximum deviation = %g at %d propensity %g \n",
   	    max_deviation,max_i,propensity[max_i]);
}
/* ---------------------------------------------------------------------- */

void AppTest::set_event(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal event command");

  fprintf(screen,
	  "Number of events: %s. Number of dependencies: %s. \n",
	  arg[0], arg[1]);
  
  n_event_types ++;
  nevents = atoi(arg[0]);
  ndep = atoi(arg[1]);
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

   int nmax = ndep;

   for (m = 0; m < nevents; m++) 
     ndepends[m] = static_cast<int>((ndep+1)*random->uniform());

   for (m = 0; m < nevents; m++)
     nmax = MAX(nmax,ndepends[m]);

   nmax = MAX(nmax,1);

   depends = memory->create_2d_int_array(nevents,nmax,
					      "test:depends");
   // select the dependencies

   for (m = 0; m < nevents; m++)
     for(int e = 0; e < ndepends[m]; e++) 
       depends[m][e] = static_cast<int>(nevents*random->uniform());

   // uncomment to print dependency graph
   //print_depend_graph();
}
/* ----------------------------------------------------------------------
   print dependency network
------------------------------------------------------------------------- */
void AppTest::print_depend_graph()
{

  int m;

  fprintf(screen, "Dependency graph: \n");
  for (m = 0; m < nevents; m++){
    fprintf(screen, "event %d  ",m);
    fprintf(screen, "ndepends %d: ",ndepends[m]);
    for(int e = 0; e < ndepends[m]; e++) {
      depends[m][e] = static_cast<int>(nevents*random->uniform());
      fprintf(screen, "%d  ", depends[m][e]);
    }
    fprintf(screen, "\n");
  }
}
/* ----------------------------------------------------------------------
   compute propensity of a single event
------------------------------------------------------------------------- */

double AppTest::compute_propensity(int m)
{
  //uniform
  //double p=.1;
  //random uniform
  //double p=random->uniform();
  //even/odd
    double p = .5 - .1 * static_cast<double>(m % 3);
  //linear
  //double p = 21.0-(double)(m+1);

  return p;
}
