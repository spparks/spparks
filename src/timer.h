/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_TIMER_H
#define SPK_TIMER_H

#include "pointers.h"

namespace SPPARKS_NS {

enum{TIME_LOOP,TIME_SOLVE,TIME_UPDATE,TIME_COMM,TIME_OUTPUT,TIME_APP,TIME_N};

class Timer : protected Pointers {
 public:
  double *array;

  Timer(class SPPARKS *);
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
