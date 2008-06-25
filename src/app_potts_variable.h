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
  void site_pick_random(int, double);
  void site_pick_local(int, double);
  double site_propensity(int, int);
  void site_event(int, int);
  void site_clear_mask(char *, int);

 private:
  int nspins;
  int *sites;
  int *spin;
};

}

#endif
