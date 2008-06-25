/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SHELL_H
#define SHELL_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Shell : protected SysPtr {
 public:
  explicit Shell(class SPPARKS *);
  void command(int, char **);
};

}

#endif
