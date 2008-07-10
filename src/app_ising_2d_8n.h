/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_ISING_2D_8N_H
#define APP_ISING_2D_8N_H

#include "app_lattice2d.h"

namespace SPPARKS_NS {

class AppIsing2d8n : public AppLattice2d {
 public:
  AppIsing2d8n(class SPPARKS *, int, char **);
  ~AppIsing2d8n();

  double site_energy(int, int);
  void site_event_rejection(int, int, class RandomPark *);
  double site_propensity(int, int);
  void site_event(int, int, int, class RandomPark *);
};

}

#endif
