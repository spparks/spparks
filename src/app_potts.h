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
AppStyle(potts,AppPotts)

#else

#ifndef SPK_APP_POTTS_H
#define SPK_APP_POTTS_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppPotts : public AppLattice {
 public:
  AppPotts(class SPPARKS *, int, char **);
  virtual ~AppPotts();

  double site_energy(int);
  virtual void site_event_rejection(int, class RandomPark *);
  virtual double site_propensity(int);
  virtual void site_event(int, class RandomPark *);

 protected:
  int nspins;
  int *sites,*unique;
};

}

#endif
#endif
