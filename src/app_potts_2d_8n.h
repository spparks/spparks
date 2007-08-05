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
  double site_propensity(int, int, int);
  void site_event(int, int, int);
  void site_update_ghosts(int, int);
  void site_clear_mask(char **, int, int);
  void survey_neighbor(const int&, const int&, int&, int[], int[]) const;

 private:
  int nspins;
};

}

#endif
