/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_NEXT_EVENT_TREE_SEARCH_H
#define SOLVE_NEXT_EVENT_TREE_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveNextEventTreeSearch : public Solve {
 public:
  SolveNextEventTreeSearch(class SPK *, int, char **);
  ~SolveNextEventTreeSearch();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  int event(double *);
  void sum_tree();
  void free_arrays();
  void set(int , double);
  int find(double);

 private:
  class RandomPark *random;
  int nevents;
  double sum;

  double *tree;
  int offset;

  int allocated;
};

}

#endif
