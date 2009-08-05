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

#ifndef SPK_SOLVE_LINEAR_H
#define SPK_SOLVE_LINEAR_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveLinear : public Solve {
 public:
  SolveLinear(class SPPARKS *, int, char **);
  ~SolveLinear();
  SolveLinear *clone();

  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  int nevents;
  double *prob;
};

}

#endif
