/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_3D_26N_H
#define APP_POTTS_3D_26N_H

#include "app_potts_3d.h"

namespace SPPARKS {

class AppPotts3d26n : public AppPotts3d {
 public:
  AppPotts3d26n(class SPK *, int, char **);
  ~AppPotts3d26n();

  double site_energy(int, int, int);
  int site_pick_random(int, int, int, double);
  int site_pick_local(int, int, int, double);
  double site_propensity(int, int, int, int);
  void site_event(int, int, int, int);
  void site_update_ghosts(int, int, int);
  void site_clear_mask(char ***, int, int, int);
  void survey_neighbor(const int&, const int&, int&, int[], int[]) const;

  void push_connected_neighbors(int, int, int, int***, int, std::stack<int>*);
  void connected_ghosts(int, int, int, int***, Cluster*, int);
};

}

#endif
