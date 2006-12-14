/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef TIMER_H
#define TIMER_H

#include "sysptr.h"

namespace SPPARKS {

#define TIME_TOTAL       0
#define TIME_LOOP        1
#define TIME_SOLVE       2
#define TIME_COMM        3
#define TIME_OUTPUT      4

#define TIME_N           5

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

 private:
  Timer(); // Not a sane operation.
  Timer(const Timer&); // Not a sane operation.
  Timer& operator=(const Timer&); // Not a sane operation.
};

}

#endif
