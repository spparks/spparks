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

#ifdef PAIR_CLASS
PairStyle(lj/cut,PairLJCut)

#else

#include "pair.h"

namespace SPPARKS_NS {

class PairLJCut : public Pair {
 public:
  PairLJCut(class SPPARKS *);
  ~PairLJCut();
  double energy(int, int, int *, double **, int *);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

  void allocate();
};

}

#endif
