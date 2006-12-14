/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_3D_STRICT_H
#define APP_GRAIN_3D_STRICT_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrain3DStrict : public AppGrain3D {
 public:
  AppGrain3DStrict(class SPK *, int, char **);
  virtual ~AppGrain3DStrict();
  virtual void run(int, char **);

 private:
  class RandomPark ***ranlat;        // array of random number generator pointers
  // These functions override definitions in base class
  // No need to make them virtual as they are called 
  // via derived class version of virtual function run().
  void iterate();
  void sweep(int,int);
};

}

#endif
