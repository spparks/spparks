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
  void site_event_rejection(int, int, class RandomPark *);
  double site_propensity(int, int);
  void site_event(int, int, int, class RandomPark *);

 private:
  double w01,w11,mu;
  double interact[4][4];

  void input_app(char *, int, char **);
};

}

#endif
