/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include <iostream>
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
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppTest::AppTest(SPK *spk, int narg, char **arg) : App(spk, narg, arg)
{
  if (narg != 1) error->all("Invalid app_style test command");

  propensity = NULL;
  ndepends = NULL;
  old_p = NULL;
  count = NULL;
  depends = NULL;
  ran_dep = NULL;

  nevents = 0;
  n_event_types = 0;
  dep_graph_flag = true;
  ntimestep = 0;
  time = 0.0;
  stoptime = 0.0;
  stats_delta = 0.0;

  // classes needed by this app
  int seed = 123124;
  random = new RandomPark(seed);

  delete timer;
  timer = new Timer(spk);
}

/* ---------------------------------------------------------------------- */

AppTest::~AppTest()
{
  delete [] propensity;
  delete [] ndepends;
  delete [] old_p;
  delete [] count;
  memory->destroy_2d_int_array(depends);
  delete [] ran_dep;

  delete random;
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

  if(dep_graph_flag){
    ndepends = new int[nevents];
    build_dependency_graph();
  }
  else{
    ran_dep = new int[2*ndep];
  }
//   if(old_p != NULL)  delete [] old_p;
//   old_p = new double[ndep];

  // compute initial propensity for each event
  // inform Nfold solver

  psum = 0;
  if (propensity != NULL) delete [] propensity;
  propensity = new double[nevents];
  for (int m = 0; m < nevents; m++) {
    propensity[m] = compute_propensity(m);
    psum += propensity[m];
  }

  //  solve->init(nevents,propensity);
  
  // allocate and zero event stats
  if (count != NULL) delete [] count;
  count = new int[2*nevents];
  for (int t = 0; t < 2*nevents; t++) count[t] = 0;

  // print stats header

  if (screen) {
    fprintf(screen,"Step Time Counts");
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"Step Time Counts");
    fprintf(logfile,"\n");
  }
  //stats();

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
  int rdp; // random number of dependencies
  int tdep; 

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {

    // uncomment to control total number of events
    nev++;

    ntimestep++;
    timer->stamp();
    ievent = solve->event(&dt);
    count[ievent]++;
    timer->stamp(TIME_SOLVE);

    // update propensity table
    // inform solver of changes

    //update event propensity

    propensity[ievent] = compute_propensity(ievent);
    solve->update(ievent, propensity);

    //update dependencies
    if (dep_graph_flag){
      for (d = 0; d < ndepends[ievent]; d++)
	propensity[depends[ievent][d]] = compute_propensity(depends[ievent][d]);
      //solve->update(ndepends[ievent],depends[ievent],old_p);
      solve->update(ndepends[ievent],depends[ievent],propensity);
    }
    else {
      rdp = 0.5 * ndep + ndep * random->uniform(); // random number of deps
      for (d = 0; d < rdp; d++){                   // for each dependency
	tdep = nevents * random->uniform();        // pick it
        ran_dep[d] = tdep;                         // record it in the array
	propensity[tdep] = compute_propensity(tdep); // calculate it
      }
      solve->update(rdp,ran_dep,propensity);       // update solver
    }

    timer->stamp(TIME_UPDATE);
    // update time by dt

    time += dt;
    if (time >= stoptime) done = 1;
    else if (ievent < 0) done = 1;
    // uncomment to control total number of events
    else if (nev > 100000) done = 1;

    if (time > stats_time || done) {
      timer->stamp();
      // comment out for huge number of events
      //stats();
      stats_time += stats_delta;
      timer->stamp(TIME_OUTPUT);
    }
  }

  timer->barrier_stop(TIME_LOOP);

  // comment out when done with timings
  printf("\nNumber of reactions = %d\n",nevents);
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
  if (ssum == 0) ssum = 1;

  if (screen) {
    for (i = 0; i < nevents; i++){
      deviation = ((double)count[i]/ssum - propensity[i]/psum)/
	((double)count[i]/ssum + propensity[i]/psum);
      if(deviation > max_deviation){
	max_deviation = deviation;
	max_i = i;
      }
    }
  }
  printf("PSUM %g %d\n",psum,ssum);
  if (screen) {
    fprintf(screen,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(screen,"%g ",(double)count[i]/ssum - propensity[i]/psum);
    fprintf(screen,"\n");
    
  }
  if (logfile) {
    fprintf(logfile,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(logfile,"%g ",(double)count[i]/ssum - propensity[i]/psum);
    fprintf(logfile,"\n");
  }

  if(max_i>nevents-1) max_i = nevents - 1;
  //  if (screen)
  // fprintf(screen,"Maximum deviation = %g at %d propensity %g \n",
  //	    max_deviation,max_i,propensity[max_i]);
}
/* ---------------------------------------------------------------------- */

void AppTest::set_event(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal event command");

  if (narg > 3)
    if (strcmp(arg[2],"lo_mem")==0){
      dep_graph_flag = false;
      random->init(atoi(arg[3]));
    }
    else error->all("Illegal event command");

  n_event_types++;
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
    for(int e = 0; e < ndepends[m]; e++){
      depends[m][e] = m;
      while (depends[m][e] == m)
	depends[m][e] = static_cast<int>(nevents*random->uniform());
    }
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
  double p = random->uniform();
  //even/odd
  //double p = .5 - .1 * static_cast<double>(m % 2);
  //linear
  //double p = (double)(m+1)/10.0;

  return p;
}

