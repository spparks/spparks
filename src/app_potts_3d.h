/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_POTTS_3D_H
#define APP_POTTS_3D_H

#include "app_lattice3d.h"

namespace SPPARKS_NS {

class AppPotts3d : public AppLattice3d {
 public:
  AppPotts3d(class SPPARKS *, int&, char **);
  ~AppPotts3d();

 protected:
  int nspins;
  char * spinfile;
};

}

#endif
