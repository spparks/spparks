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
AppStyle(potts/ca,AppPottsCA)

#else

#ifndef SPK_APP_POTTS_CA_H
#define SPK_APP_POTTS_CA_H

#include "app_potts_neighonly.h"

namespace SPPARKS_NS {

class AppPottsCA : public AppPottsNeighOnly {
 public:
  AppPottsCA(class SPPARKS *, int, char **);
  ~AppPottsCA();
  void grow_app();
  void init_app();
  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  void app_update(double);
  //void iterate_rejection(double);
  
 private:
  double *stored_e;
  int *nuc_switch;
  
  double rc_repeats;
  double gg_repeats;
  double nuclei_pct;
  double threshold;
  double chance;
  
  int aspect;
  int counter;

  class RandomPark *ranlatt;
};
	
}

#endif
#endif
