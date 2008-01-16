/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_ISING_H
#define APP_ISING_H

#include "app_lattice.h"

namespace SPPARKS {

class AppIsing : public AppLattice {
 public:
  AppIsing(class SPK *, int, char **);
  ~AppIsing();

  double site_energy(int);
  void site_pick_random(int, double);
  void site_pick_local(int, double);
  double site_propensity(int, int);
  void site_event(int, int);
  void site_clear_mask(char *, int);

 private:
  int *sites;
};

}

#endif
