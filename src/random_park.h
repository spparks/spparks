/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef RANDOM_PARK_H
#define RANDOM_PARK_H

namespace SPPARKS {

class RandomPark {
 public:
  int seed;

  RandomPark(int);
  ~RandomPark() {}
  void init(int);
  double uniform();
  int irandom(int);
};

}

#endif
