/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_TREE_H
#define SOLVE_TREE_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveTree : public Solve {
 public:
  SolveTree(class SPPARKS *, int, char **);
  ~SolveTree();
  SolveTree *clone();

  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);
  void sum_tree();
  void free_arrays();
  void set(int, double);
  int find(double);
  void tree_to_screen(int);

 private:
  int seed;
  class RandomPark *random;
  int nevents;
  int ntotal;

  double *tree;
  int offset;
  int allocated;
};

}

#endif
