/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef ERROR_H
#define ERROR_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Error : protected SysPtr {
 public:
  explicit Error(class SPPARKS *);

  void universe_all(const char *);
  void universe_one(const char *);

  void all(const char *);
  void one(const char *);
  void warning(const char *);
};

}

#endif
