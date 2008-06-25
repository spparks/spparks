/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef RANDOM_FIBLAG_H
#define RANDOM_FIBLAG_H

namespace SPPARKS_NS {

class RandomFibLag {
 public:
  RandomFibLag(int);
  ~RandomFibLag() {}
  double uniform();
  int irandom(int);

 private:
  int seed;
  double fib[56];
  int i0,i8,i16,i24,i55;
};

}

#endif
