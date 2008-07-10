/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_VARIABLE_H
#define APP_POTTS_VARIABLE_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppPottsVariable : public AppLattice {
 public:
  AppPottsVariable(class SPPARKS *, int, char **);
  ~AppPottsVariable();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int nspins;
  int *sites;
  int *spin;
};

}

#endif
