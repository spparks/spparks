/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef FINISH_H
#define FINISH_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Finish : protected SysPtr {
 public:
  explicit Finish(class SPPARKS *);
  ~Finish() {}

 private:
  void stats(int, double *, double *, double *, double *, int, int *);
};

}

#endif
