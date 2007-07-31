/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_ISING_3D_6N_H
#define APP_ISING_3D_6N_H

#include "app_lattice3d.h"

namespace SPPARKS {

class AppIsing3d6n : public AppLattice3d {
 public:
  AppIsing3d6n(class SPK *, int, char **);
  ~AppIsing3d6n();

  double site_energy(int, int, int);
  int site_pick_random(int, int, int, double);
  int site_pick_local(int, int, int, double);
  double site_propensity(int, int, int, int);
  void site_event(int, int, int, int);
  void site_update_ghosts(int, int, int);
  void site_clear_mask(char ***, int, int, int);
};

}

#endif
