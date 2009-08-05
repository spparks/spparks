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

#ifndef SPK_SOLVE_H
#define SPK_SOLVE_H

#include "pointers.h"

namespace SPPARKS_NS {

class Solve : protected Pointers {
 public:
  char *style;

  Solve(class SPPARKS *, int, char **);
  virtual ~Solve();

  double get_total_propensity();
  int get_num_active();

  // pure virtual functions, must be defined in child class

  virtual Solve *clone() = 0;

  virtual void init(int, double *) = 0;
  virtual void update(int, int *, double *) = 0;
  virtual void update(int, double *) = 0;
  virtual void resize(int, double *) = 0;
  virtual int event(double *) = 0;

 protected:
  double sum;
  int num_active;
};

}

#endif
