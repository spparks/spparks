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

#ifndef SPK_RANDOM_MARS_H
#define SPK_RANDOM_MARS_H

#include "pointers.h"

namespace SPPARKS_NS {

class RanMars : protected Pointers {
 public:
  RanMars(class SPPARKS *);
  ~RanMars();
  void init(int);
  double uniform();

 private:
  int initflag;
  int i97,j97;
  double c,cd,cm;
  double *u;
};

}

#endif

/* ERROR/WARNING messages:

E: Seed command has not been used

The seed command must be used if another command requires random
numbers.

*/
