/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_2D_24N_H
#define APP_POTTS_2D_24N_H

#include "app_potts_2d.h"

namespace SPPARKS_NS {

class AppPotts2d24n : public AppPotts2d {
 public:
  AppPotts2d24n(class SPPARKS *, int, char **);
  ~AppPotts2d24n();

  double site_energy(int, int);
  void site_event_rejection(int, int, class RandomPark *);
  double site_propensity(int, int);
  void site_event(int, int, int, class RandomPark *);
};

}

#endif
