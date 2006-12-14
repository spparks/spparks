/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_SWEEP_H
#define SOLVE_SWEEP_H

#include "solve.h"

namespace SPPARKS {

class SolveSweep : public Solve {
 public:
  SolveSweep(class SPK *, int, char **);
  ~SolveSweep();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  int event(double *);
};

}

#endif
