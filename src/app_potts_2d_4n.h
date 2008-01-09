/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_2D_4N_H
#define APP_POTTS_2D_4N_H

#include "app_potts_2d.h"

namespace SPPARKS {

class AppPotts2d4n : public AppPotts2d {
 public:
  AppPotts2d4n(class SPK *, int, char **);
  ~AppPotts2d4n();

  double site_energy(int, int);
  int site_pick_random(int, int, double);
  int site_pick_local(int, int, double);
  double site_propensity(int, int, int);
  void site_event(int, int, int);
  void site_update_ghosts(int, int);
  void site_clear_mask(char **, int, int);
  void survey_neighbor(const int&, const int&, int&, int[], int[]) const;

};

}

#endif
