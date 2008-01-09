/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_2D_H
#define APP_POTTS_2D_H

#include "app_lattice2d.h"

namespace SPPARKS {

class AppPotts2d : public AppLattice2d {
 public:
  AppPotts2d(class SPK *, int&, char **);
  ~AppPotts2d();

 protected:
  int nspins;
  char * spinfile;
};

}

#endif
