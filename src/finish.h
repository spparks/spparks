/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef FINISH_H
#define FINISH_H

#include "sysptr.h"

namespace SPPARKS {

class Finish : protected SysPtr {
 public:
  explicit Finish(class SPK *);
  ~Finish() {}

 private:
  void stats(int, double *, double *, double *, double *, int, int *);

 private:
  Finish(); // Not a sane operation.
  Finish(const Finish&); // Not a sane operation.
  Finish& operator=(const Finish&); // Not a sane operation.
};

}

#endif
