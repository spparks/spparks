/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_MEMBRANE_H
#define APP_MEMBRANE_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppMembrane : public AppLattice {
 public:
  AppMembrane(class SPPARKS *, int, char **);
  ~AppMembrane();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  double w01,w11,mu;
  double interact[4][4];
  int *sites;

  void input_app(char *, int, char **);
};

}

#endif
