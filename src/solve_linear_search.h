/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_LINEAR_SEARCH_H
#define SOLVE_LINEAR_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveLinearSearch : public Solve {
 public:
  SolveLinearSearch(class SPK *, int, char **);
  ~SolveLinearSearch();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  int nevents,nzeroes;
  double *prob;
  double sum;
};

}

#endif
