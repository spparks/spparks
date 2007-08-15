/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_LINEAR_H
#define SOLVE_LINEAR_H

#include "solve.h"

namespace SPPARKS {

class SolveLinear : public Solve {
 public:
  SolveLinear(class SPK *, int, char **);
  ~SolveLinear();
  SolveLinear *clone();

  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  int seed;
  class RandomPark *random;
  int nevents,nzeroes;
  double *prob;
  double sum;
};

}

#endif
