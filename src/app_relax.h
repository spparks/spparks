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

#ifdef APP_CLASS
AppStyle(relax,AppRelax)

#else

#include "app_off_lattice.h"

namespace SPPARKS_NS {

class AppRelax : public AppOffLattice {
 public:
  AppRelax(class SPPARKS *, int, char **);
  ~AppRelax();
  void grow_app();
  void init_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int) {return 0.0;}
  void site_event(int, class RandomPark *) {}

 private:
  double delta,deltasq;
  int *site;
  class Pair *pair;
  
  double site_energy_neighbor(int);
};

}

#endif
