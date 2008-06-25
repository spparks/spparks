/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_MEMBRANE_2D_H
#define APP_MEMBRANE_2D_H

#include "app_lattice2d.h"

namespace SPPARKS_NS {

class AppMembrane2d : public AppLattice2d {
 public:
  AppMembrane2d(class SPPARKS *, int, char **);
  ~AppMembrane2d();

  double site_energy(int, int);
  int site_pick_random(int, int, double);
  int site_pick_local(int, int, double);
  double site_propensity(int, int, int);
  void site_event(int, int, int);
  void site_update_ghosts(int, int);
  void site_clear_mask(char **, int, int);

 private:
  double w01,w11,mu;
  double interact[4][4];

  void input_app(char *, int, char **);
};

}

#endif
