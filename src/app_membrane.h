/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_MEMBRANE_H
#define APP_MEMBRANE_H

#include "app_lattice.h"

namespace SPPARKS {

class AppMembrane : public AppLattice {
 public:
  AppMembrane(class SPK *, int, char **);
  ~AppMembrane();

  double site_energy(int);
  int site_pick_random(int, double);
  int site_pick_local(int, double);
  double site_propensity(int, int);
  void site_event(int, int);
  void site_clear_mask(char *, int);

 private:
  double w01,w11,mu;
  double interact[4][4];
  int *sites;

  void input_app(char *, int, char **);
};

}

#endif
