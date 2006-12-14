/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SHELL_H
#define SHELL_H

#include "sysptr.h"

namespace SPPARKS {

class Shell : protected SysPtr {
 public:
  explicit Shell(class SPK *);
  void command(int, char **);

 private:
  Shell(); // Not a sane operation.
  Shell(const Shell&); // Not a sane operation.
  Shell& operator=(const Shell&); // Not a sane operation.
};

}

#endif
