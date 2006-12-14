/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_TEST_H
#define APP_TEST_H

#include "app.h"

namespace SPPARKS {

class AppTest : public App {
 public:
  AppTest(class SPK *, int, char **);
  // virtual keyword informative but not required
  // since these functions are virtual in base class
  virtual ~AppTest();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);

 private:

  //time keeping
  int ntimestep;
  double time,stoptime;

  //event table
  int *ndepends;             // # of events that depend on each event
  int **depends;             // i,j = jth event that depends on ith event
  double *propensity;        // propensity of each event
  int nevents;
  int n_event_types;

  //event properties

  //stats timing
  double stats_time,stats_delta;

  //stats
  int count;

  //methods
  void iterate();
  void stats();

  void set_stats(int, char **);

  void build_dependency_graph();
  double compute_propensity(int);

  void set_event(int, char **);
};

}

#endif
