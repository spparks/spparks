/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GROUP2_H
#define SOLVE_GROUP2_H

#include "solve.h"

namespace SPPARKS {

class SolveGroup2 : public Solve {
 public:
  SolveGroup2(class SPK *, int, char **);
  ~SolveGroup2();
  SolveGroup2 *clone();

  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  int seed;
  class RandomPark *random;
  class Groups2 *groups;
  int nevents,nzeroes;

  double *p;                     // local copy of propensities
  double sum;
  int ngroups_in;
  bool ngroups_flag;
  double lo,hi;
};

}

#endif

