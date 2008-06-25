/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GPPT
#define SOLVE_GPPT

#include "solve.h"
#include "node.h"

namespace SPPARKS_NS {
  
class SolveGPPT : public Solve {
 public:
  SolveGPPT(class SPPARKS *, int, char **);
  ~SolveGPPT();
  SolveGPPT *clone();

  void input(int, char **) {}
  void init(int, double *) {}
  void update(int, double *){}
  void update(int, int *, double *){}
  void resize(int, double *);
  int event(double *);
    
 private:
  int seed;
  class RandomPark *random;
};
  
}

#endif
