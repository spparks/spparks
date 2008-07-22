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

#ifndef RANDOM_RANECU_H
#define RANDOM_RANECU_H

namespace SPPARKS_NS {

class RandomRanecu {
 public:
  RandomRanecu(int, int);
  ~RandomRanecu() {}
  double uniform();
  int irandom(int);

 private:
  int seed1,seed2;
};

}

#endif
