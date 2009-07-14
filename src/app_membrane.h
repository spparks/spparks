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

#ifndef APP_MEMBRANE_H
#define APP_MEMBRANE_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppMembrane : public AppLattice {
 public:
  AppMembrane(class SPPARKS *, int, char **);
  ~AppMembrane();
  void input_app(char *, int, char **);

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  double w01,w11,mu;
  double interact[4][4];
  int *sites;
};

}

#endif
