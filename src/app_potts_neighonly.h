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

#ifndef APP_POTTS_NEIGHONLY_H
#define APP_POTTS_NEIGHONLY_H

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsNeighOnly : public AppPotts {
 public:
  AppPottsNeighOnly(class SPPARKS *, int, char **);
  ~AppPottsNeighOnly() {}

  void site_event_rejection(int, class RandomPark *);
};

}

#endif
