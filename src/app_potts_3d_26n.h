/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_3D_26N_H
#define APP_POTTS_3D_26N_H

#include "app_potts_3d.h"

namespace SPPARKS_NS {

class AppPotts3d26n : public AppPotts3d {
 public:
  AppPotts3d26n(class SPPARKS *, int, char **);
  ~AppPotts3d26n();

  double site_energy(int, int, int);
  void site_event_rejection(int, int, int, class RandomPark *);
  double site_propensity(int, int, int);
  void site_event(int, int, int, int, class RandomPark *);

  void push_connected_neighbors(int, int, int, int***, int, std::stack<int>*);
  void connected_ghosts(int, int, int, int***, Cluster*, int);
};

}

#endif
