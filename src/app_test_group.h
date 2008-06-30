/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_TEST_GROUP_H
#define APP_TEST_GROUP_H

#include "app.h"

namespace SPPARKS_NS {

class AppTestGroup : public App {
 public:
  AppTestGroup(class SPPARKS *, int, char **);
  ~AppTestGroup();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

 private:
  class RandomPark *random;

  int ntimestep;
  double time;
  int nlimit;

  int nevents;               // # of events
  double *propensity;        // propensity of each event
  double tweak;              // size of propensity tweak
  double psum;
  int *count;

  int dep_graph_flag;
  int ndep;                  // max number of dependencies from user
  int *ndepends;             // # of events that depend on each event
  int **depends;             // i,j = jth event that depends on ith event
  int *ran_dep;              // random deps for on-the-fly generation

  void iterate();
  void build_dependency_graph();
  double compute_propensity(int);
  void set_event(int, char **);

  void stats(char *);
  void stats_header(char *);
};

}

#endif
