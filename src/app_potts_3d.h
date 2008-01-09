/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_3D_H
#define APP_POTTS_3D_H

#include "app_lattice3d.h"

namespace SPPARKS {

class AppPotts3d : public AppLattice3d {
 public:
  AppPotts3d(class SPK *, int&, char **);
  ~AppPotts3d();

 protected:
  int nspins;
  char * spinfile;
};

}

#endif
