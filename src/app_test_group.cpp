/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include "string.h"
#include "app_test_group.h"
#include "spk.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "math.h"

using namespace SPPARKS;
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppTestGroup::AppTestGroup(SPK *spk, int narg, char **arg) :
  App(spk, narg, arg)
{
  if (narg != 1) error->all("Invalid app_style test/group command");

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
  stats_delta = 0.0;

  // classes needed by this app
  int seed = 123124;
  random = new RandomPark(seed);

  delete timer;
  timer = new Timer(spk);
}

/* ---------------------------------------------------------------------- */

AppTestGroup::~AppTestGroup()
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

void AppTestGroup::init()
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
    double tp = pow(2,-random->uniform()*20);
    //double tp = random->uniform();
    if(tp > 1.0e-6) propensity[m] = tp;
    else propensity[m] = 1.0e-6;
    psum += propensity[m];
  }

  //  solve->init(nevents,propensity);
  
  // allocate and zero event stats
  if (count != NULL) delete [] count;
  count = new int[nevents];
  for (int t = 0; t < nevents; t++) count[t] = 0;

  // print stats header

  if (screen) {
    fprintf(screen,"Step Time Counts");
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf(logfile,"Step Time Counts");
    fprintf(logfile,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void AppTestGroup::input(char *command, int narg, char **arg)
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

void AppTestGroup::run(int narg, char **arg)
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

void AppTestGroup::iterate()
{
  int d,ievent;
  double dt;
  int rdp; // random number of dependencies
  int tdep; 

  ntimestep = 0;
  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
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
    if (ntimestep >= nlimit) done = 1;
    else if (ievent < 0) done = 1;
  }

  timer->barrier_stop(TIME_LOOP);
  
  printf("\nNumber of reactions, events = %d %d\n",nevents,ntimestep);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppTestGroup::stats()
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
      fprintf(screen,"%6.3d ",(double)count[i]/ssum - propensity[i]/psum);
    fprintf(screen,"\n");
    
  }
  if (logfile) {
    fprintf(logfile,"%d %g ",ntimestep,time);
    for (i = 0; i < nevents; i++)
      fprintf(logfile,"%6.3d ",(double)count[i]/ssum - propensity[i]/psum);
    fprintf(logfile,"\n");
  }

  if(max_i>nevents-1) max_i = nevents - 1;
  //  if (screen)
  // fprintf(screen,"Maximum deviation = %g at %d propensity %g \n",
  //	    max_deviation,max_i,propensity[max_i]);
}
/* ---------------------------------------------------------------------- */

void AppTestGroup::set_event(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal event command");

  if (narg > 4)
    if (strcmp(arg[3],"lo_mem")==0){
      dep_graph_flag = false;
      random->init(atoi(arg[4]));
    }
    else error->all("Illegal event command");

  n_event_types++;
  nevents = atoi(arg[0]);
  ndep = atoi(arg[1]);
  tweak = atof(arg[2]); if(tweak>0) tweak /= 100.0;
}

/* ---------------------------------------------------------------------- */

void AppTestGroup::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
}

/* ----------------------------------------------------------------------
   build dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void AppTestGroup::build_dependency_graph()
{
  int m;
  
  int nmax = ndep;
  
  for (int m = 0; m < nevents; m++) 
    ndepends[m] = static_cast<int>(ndep*random->uniform()) + 1;
  //for (m = 0; m < nevents; m++) 
  //  ndepends[m] = static_cast<int>((ndep+1)*random->uniform());

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

void AppTestGroup::print_depend_graph()
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

double AppTestGroup::compute_propensity(int m)
{
  double p;
  //uniform
  //double p=.1;
  //random uniform tweak
  p = propensity[m];
  p += p*tweak*(random->uniform()-0.5);
  
  //even/odd
  //double p = .5 - .1 * static_cast<double>(m % 2);
  //linear
  //double p = (double)(m+1)/10.0;
  if(p>1.0) p = 0.99999;
  else if(p<1e-6) p = 1e-6;

  return p;
}

