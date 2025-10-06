/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_FINISH_H
#define SPK_FINISH_H

#include "pointers.h"

namespace SPPARKS_NS {

class Finish : protected Pointers {
 public:
  Finish(class SPPARKS *, int);
  ~Finish() {}

 private:
  void stats(int, double *, double *, double *, double *, int, int *);
};

}

#endif
