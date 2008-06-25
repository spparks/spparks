/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_H
#define DIAG_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Diag : protected SysPtr {
 public:
  Diag(class SPPARKS *, int, char **);
  virtual ~Diag();

  int check_time(double, int);
  void setup_time(double);

  // virtual functions with default definitions
  // may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};

  // pure virtual functions, must be defined in child class
  
  virtual void init(double) = 0;
  virtual void compute(double, int) = 0;

 protected:
  char *style;
  int me,nprocs;

  double diag_time,diag_delta,diag_scale,diag_t0;
  int diag_nrepeat,diag_irepeat,diag_ilogfreq;

};

}

#endif
