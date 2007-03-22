/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef TIMER_H
#define TIMER_H

#include "sysptr.h"

namespace SPPARKS {

  enum{TIME_LOOP,TIME_SOLVE,TIME_UPDATE,TIME_COMM,TIME_OUTPUT,TIME_APP,TIME_N};

class Timer : protected SysPtr {
 public:
  double *array;

  explicit Timer(class SPK *);
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

 private:
  double previous_time;
};

}

#endif
