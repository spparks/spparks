/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GROUP_H
#define SOLVE_GROUP_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveGroup : public Solve {
 public:
  SolveGroup(class SPPARKS *, int, char **);
  ~SolveGroup();
  SolveGroup *clone();

  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  int seed;
  class RandomPark *random;
  class Groups *groups;
  int nevents;

  double *p;                     // local copy of propensities
  double sum;
  double lo,hi;
  int ngroups;
};

}

#endif

