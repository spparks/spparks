/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_H
#define SOLVE_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Solve : protected SysPtr {
 protected:
  double sum;
  int num_active;

 public:
  char *style;

  Solve(class SPPARKS *, int, char **);
  virtual ~Solve();

  double get_total_propensity();
  int get_num_active();

  // pure virtual functions, must be defined in child class

  virtual Solve *clone() = 0;

  virtual void input(int, char **) = 0;
  virtual void init(int, double *) = 0;
  virtual void update(int, int *, double *) = 0;
  virtual void update(int, double *) = 0;
  virtual void resize(int, double *) = 0;
  virtual int event(double *) = 0;
};

}

#endif
