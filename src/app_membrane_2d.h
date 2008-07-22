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

#ifndef APP_MEMBRANE_2D_H
#define APP_MEMBRANE_2D_H

#include "app_lattice2d.h"

namespace SPPARKS_NS {

class AppMembrane2d : public AppLattice2d {
 public:
  AppMembrane2d(class SPPARKS *, int, char **);
  ~AppMembrane2d();

  double site_energy(int, int);
  void site_event_rejection(int, int, class RandomPark *);
  double site_propensity(int, int);
  void site_event(int, int, int, class RandomPark *);

 private:
  double w01,w11,mu;
  double interact[4][4];

  void input_app(char *, int, char **);
};

}

#endif
