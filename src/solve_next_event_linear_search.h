/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_NEXT_EVENT_LINEAR_SEARCH_H
#define SOLVE_NEXT_EVENT_LINEAR_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveNextEventLinearSearch : public Solve {
 public:
  SolveNextEventLinearSearch(class SPK *, int, char **);
  ~SolveNextEventLinearSearch();
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
