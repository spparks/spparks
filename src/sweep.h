/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_H
#define SWEEP_H

#include "sysptr.h"

namespace SPPARKS {

class Sweep : protected SysPtr {
 public:
  char *style;

  Sweep(class SPK *, int, char **);
  virtual ~Sweep();

  // pure virtual functions, must be defined in child class

  virtual void input(int, char **) = 0; 
  virtual void do_sweep(double&) = 0;
  virtual double compute_energy() = 0;
};

}

#endif
