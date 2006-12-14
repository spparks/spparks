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

  void universe_all(char *);
  void universe_one(char *);

  void all(char *);
  void one(char *);
  void warning(char *);

 private:
  Error(); // Not a sane operation.
  Error(const Error&); // Not a sane operation.
  Error& operator=(const Error&); // Not a sane operation.
};

}

#endif
