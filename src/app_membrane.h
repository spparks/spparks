/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_MEMBRANE_H
#define APP_MEMBRANE_H

#include "app_lattice2d.h"

namespace SPPARKS {

class AppMembrane : public AppLattice2d {
 public:
  AppMembrane(class SPK *, int, char **);
  ~AppMembrane();

  double site_energy(int, int);
  int site_pick_random(int, int, double);
  int site_pick_local(int, int, double);
  double site_propensity(int, int, int);
  void site_event(int, int, int);
  void site_update_ghosts(int, int);
  void site_clear_mask(char **, int, int);

 private:
  double w01,w11,w22,prefactor;

  void input_app(char *, int, char **);
};

}

#endif
