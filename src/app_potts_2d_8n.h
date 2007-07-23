/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_2D_8N_H
#define APP_POTTS_2D_8N_H

#include "app_lattice2d.h"

namespace SPPARKS {

class AppPotts2d8n : public AppLattice2d {
 public:
  AppPotts2d8n(class SPK *, int, char **);
  ~AppPotts2d8n();

  double site_energy(int, int);
  int site_pick_random(int, int, double);
  int site_pick_local(int, int, double);
  double site_propensity(int, int);
  void site_event(int, int);
  void site_event_sector(int, int);
  void site_update_ghost(int, int);
  void site_clear_mask(char **, int, int);

 private:
  int nspins;
};

}

#endif
