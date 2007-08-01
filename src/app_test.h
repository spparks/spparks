/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

/*----------------------------------------------------------------------
This class tests the various ways of generating the stochastic time series
from a non-stationary discrete event probability distribution.

Blocks:
0. Read user input
1. Generate initial discrete distribution by compiling a table of possible
events and their respective propensities.
2. Pick an event to happen next.
3. Update distribution.
4. Increment time.
5. Print stats(?).
6. Go to 2.

Distribution properties to test:

0. Propensity table size

1. Propensity table generation:
       a) Random uniform
       b) Random non-uniform

2. Propensity table updates:
       a) Independent
       b) Dependent with variable number of dependencies

Questions:

1. When does an update cost roughly the same as an initialization?

Code changes:

1. input.cpp
added   else if (!strcmp(command,"event")) event();
1. style_user.h
added definitions for the test app and next_event_linear_search style
----------------------------------------------------------------------*/

#ifndef APP_TEST_H
#define APP_TEST_H

#include "app.h"

namespace SPPARKS {

class AppTest : public App {
 public:
  AppTest(class SPK *, int, char **);
  ~AppTest();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

 private:

  class RandomPark *random;

  //time keeping
  int ntimestep;
  double time,stoptime;

  //event table
  int *ndepends;             // # of events that depend on each event
  int **depends;             // i,j = jth event that depends on ith event
  double *propensity;        // propensity of each event
  int nevents;
  int ndep;                  // max number of dependencies from user
  double *old_p;             // old propensity values
  int n_event_types;

  //event properties

  //stats timing
  double stats_time,stats_delta;

  //stats
  int *count;
  int ssum;
  double psum;

  //methods
  void iterate();
  void stats();

  void set_stats(int, char **);

  void build_dependency_graph();
  void print_depend_graph();
  double compute_propensity(int);

  void set_event(int, char **);

};

}

#endif

