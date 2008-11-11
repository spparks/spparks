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

#ifndef APP_POTTS_3D_26N_H
#define APP_POTTS_3D_26N_H

#include "app_lattice3d.h"

namespace SPPARKS_NS {

class AppPotts3d26n : public AppLattice3d {
 public:
  AppPotts3d26n(class SPPARKS *, int, char **);
  ~AppPotts3d26n();

  double site_energy(int, int, int);
  void site_event_rejection(int, int, int, class RandomPark *);
  double site_propensity(int, int, int);
  void site_event(int, int, int, int, class RandomPark *);

 private:
  int nspins;

  void push_connected_neighbors(int, int, int, int***, int, std::stack<int>*);
  void connected_ghosts(int, int, int, int***, Cluster*, int);
};

}

#endif
