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

#ifndef RANDOM_FIBLAG_H
#define RANDOM_FIBLAG_H

namespace SPPARKS_NS {

class RandomFibLag {
 public:
  RandomFibLag(int);
  ~RandomFibLag() {}
  double uniform();
  int irandom(int);

 private:
  int seed;
  double fib[56];
  int i0,i8,i16,i24,i55;
};

}

#endif
