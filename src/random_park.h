/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

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
