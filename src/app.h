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
  enum AppClasses {GENERAL,LATTICE,LATTICE2D,LATTICE3D};
  char *style;
  int appclass;

  App(class SPK *, int, char **);
  virtual ~App();
  
  // virtual functions with empty definitions
  // may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};

  virtual void dump_header(){};
  virtual void dump(){};
  virtual void set_stats(int, char **){};
  virtual void set_dump(int, char **){};

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
