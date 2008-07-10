/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_ISING_H
#define APP_ISING_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppIsing : public AppLattice {
 public:
  AppIsing(class SPPARKS *, int, char **);
  ~AppIsing();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int *sites;
};

}

#endif
