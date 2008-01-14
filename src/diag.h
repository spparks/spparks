/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_H
#define DIAG_H

#include "sysptr.h"

namespace SPPARKS {

class Diag : protected SysPtr {
 public:
  char *style;
  
  Diag(class SPK *, int, char **);
  virtual ~Diag();
  
  // virtual functions with empty definitions
  // may be overridden in child class

  virtual void foo(int, char **){};

  // pure virtual functions, must be defined in child class
  
  virtual void init(double) = 0;
  virtual void compute(double, int) = 0;

 private:
  Diag(); // Not a sane operation.
  Diag(const Diag&); // Not a sane operation.
  Diag& operator=(const Diag&); // Not a sane operation.

 protected:
  double diag_time,diag_delta;
  int me,nprocs;

};

}

#endif
