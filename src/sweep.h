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
  virtual void init() = 0;
  virtual void do_sweep(double&) = 0;
};

}

#endif
