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

#ifdef APP_CLASS
AppStyle(test/group,AppTestGroup)

#else

#ifndef SPK_APP_TEST_GROUP_H
#define SPK_APP_TEST_GROUP_H

#include "app.h"

namespace SPPARKS_NS {

class AppTestGroup : public App {
 public:
  AppTestGroup(class SPPARKS *, int, char **);
  ~AppTestGroup();
  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();

 private:
  class RandomPark *random;

  int ncount;                // # of events performed

  int nevents;               // # of user-defined events (reactions)
  double *propensity;        // propensity of each event
  double pmax,pmin;          // maximum/minimun propensity value
  double tweak;              // percentage propensity tweak
  int seed;                  // random number seed

  double psum;
  int *count;

  int dep_graph;             // 1 if build/store dependency graph, else 0
  int ndep;                  // max number of dependencies from user
  int *ndepends;             // # of events that depend on each event
  int **depends;             // i,j = jth event that depends on ith event
  int *ran_dep;              // random deps for on-the-fly generation

  void build_dependency_graph();
  double compute_propensity(int);

  void stats(char *);
  void stats_header(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Invalid event count for app_style test/group

Number of events must be > 0.

E: Invalid probability bounds for app_style test/group

Self-explanatory.

E: Invalid probability delta for app_style test/group

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: No solver class defined

Self-explanatory.

*/
