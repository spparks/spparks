/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GILLESPIE_H
#define SOLVE_GILLESPIE_H

#include "solve.h"

namespace SPPARKS {

class SolveGillespie : public Solve {
 public:
  SolveGillespie(class SPK *, int, char **);
  ~SolveGillespie();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  int nreactions;
  double *prob;
  double sum;
};

}

#endif
