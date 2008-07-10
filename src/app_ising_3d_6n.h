/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_ISING_3D_6N_H
#define APP_ISING_3D_6N_H

#include "app_lattice3d.h"

namespace SPPARKS_NS {

class AppIsing3d6n : public AppLattice3d {
 public:
  AppIsing3d6n(class SPPARKS *, int, char **);
  ~AppIsing3d6n();

  double site_energy(int, int, int);
  void site_event_rejection(int, int, int, class RandomPark *);
  double site_propensity(int, int, int);
  void site_event(int, int, int, int, class RandomPark *);
};

}

#endif
