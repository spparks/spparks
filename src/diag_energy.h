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

#ifndef DIAG_ENERGY_H
#define DIAG_ENERGY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagEnergy : public Diag {

 public:
  DiagEnergy(class SPPARKS *, int, char **);
  virtual ~DiagEnergy();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  class AppLattice *applattice;

  int nlocal;

  double energy;

  FILE* fp;

};

}

#endif
