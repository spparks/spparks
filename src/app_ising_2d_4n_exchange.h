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

#ifndef APP_ISING_2D_4N_EXCHANGE_H
#define APP_ISING_2D_4N_EXCHANGE_H

#include "app_lattice2d.h"

namespace SPPARKS_NS {

class AppIsing2d4nExchange : public AppLattice2d {
 public:
  AppIsing2d4nExchange(class SPPARKS *, int, char **);
  ~AppIsing2d4nExchange();

  double site_energy(int, int);
  void site_event_rejection(int, int, class RandomPark *);
  double site_propensity(int, int);
  void site_event(int, int, int, class RandomPark *);
};

}

#endif
