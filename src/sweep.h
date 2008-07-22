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

#ifndef SWEEP_H
#define SWEEP_H

#include "pointers.h"

namespace SPPARKS_NS {

class Sweep : protected Pointers {
 public:
  char *style;
  bool Lmask,Lstrict,Lkmc,Ladapt;

  Sweep(class SPPARKS *, int, char **);
  virtual ~Sweep();
  virtual void init() = 0;
  virtual void do_sweep(double&) = 0;
};

}

#endif
