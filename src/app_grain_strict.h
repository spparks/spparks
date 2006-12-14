/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_STRICT_H
#define APP_GRAIN_STRICT_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrainStrict : public AppGrain {
 public:
  AppGrainStrict(class SPK *, int, char **);
  virtual ~AppGrainStrict();
  virtual void run(int, char **);

 private:
  class RandomPark **ranlat;        // array of random number generator pointers
  // These functions override definitions in base class
  // No need to make them virtual as they are called 
  // via derived class version of virtual function run().
  void iterate();
  void sweep(int,int);
};

}

#endif
