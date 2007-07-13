/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_GPPT
#define SOLVE_GPPT

#include "solve.h"
#include "node.h"

namespace SPPARKS {
  
class SolveGPPT : public Solve {
  public:
    SolveGPPT(class SPK *, int, char **);
    ~SolveGPPT();
    void input(int, char **) {}
    void init(int, double *) {}
    void update(int, double *){}
    void update(int, int *, double *){}
    void resize(int, double *);
    int event(double *);
    
  private:
    class RandomPark *random;
    
  };
  
}

#endif
