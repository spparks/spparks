/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef ERROR_H
#define ERROR_H

#include "sysptr.h"

namespace SPPARKS {

class Error : protected SysPtr {
 public:
  explicit Error(class SPK *);

  void universe_all(const char *);
  void universe_one(const char *);

  void all(const char *);
  void one(const char *);
  void warning(const char *);

 private:
  Error(); // Not a sane operation.
  Error(const Error&); // Not a sane operation.
  Error& operator=(const Error&); // Not a sane operation.
};

}

#endif
