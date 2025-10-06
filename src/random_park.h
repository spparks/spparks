/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_RANDOM_PARK_H
#define SPK_RANDOM_PARK_H

#include "spktype.h"

namespace SPPARKS_NS {

class RandomPark {
 public:
  int seed;

  RandomPark(int);
  RandomPark(double);
  ~RandomPark() {}
  void reset(double, int, int);
  void tagreset(double, tagint, int);
  double uniform();
  int irandom(int);
  tagint tagrandom(tagint);
  bigint bigrandom(bigint);
};

}

#endif
