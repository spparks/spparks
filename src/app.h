/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_H
#define APP_H

#include "sysptr.h"

namespace SPPARKS {

class App : protected SysPtr {
 public:
  char *style;
  
  App(class SPK *, int, char **);
  virtual ~App() {}
  
  // pure virtual functions, must be defined in child class
  
  virtual void input(char *, int, char **) = 0;
  virtual void init() = 0;
  virtual void run(int, char **) = 0;

 private:
  App(); // Not a sane operation.
  App(const App&); // Not a sane operation.
  App& operator=(const App&); // Not a sane operation.
};

}

#endif
