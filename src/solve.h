/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_H
#define SOLVE_H

#include "sysptr.h"

namespace SPPARKS {

class Solve : protected SysPtr {
 public:
  char *style;

  Solve(class SPK *, int, char **);
  virtual ~Solve() {}

  // pure virtual functions, must be defined in child class

  virtual void input(int, char **) = 0;
  virtual void init(int, double *) = 0;
  virtual void update(int, int *, double *) = 0;
  virtual int event(double *) = 0;

 private:
  Solve(); // Not a sane operation.
  Solve(const Solve&); // Not a sane operation.
  Solve& operator=(const Solve&); // Not a sane operation.
};

}

#endif
