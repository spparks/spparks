/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef RANDOM_RANECU_H
#define RANDOM_RANECU_H

namespace SPPARKS_NS {

class RandomRanecu {
 public:
  RandomRanecu(int, int);
  ~RandomRanecu() {}
  double uniform();
  int irandom(int);

 private:
  int seed1,seed2;
};

}

#endif
