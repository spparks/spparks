/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GROUP_SEARCH_H
#define SOLVE_GROUP_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveGroupSearch : public Solve {
 public:
  SolveGroupSearch(class SPK *, int, char **);
  ~SolveGroupSearch();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  class Groups *groups;
  int nevents;
  double *p;
  double last_size;
  double sum;
  int seed;
  int ngroups_in;
  bool ngroups_flag;
  double lo, hi;
};

}

#endif

