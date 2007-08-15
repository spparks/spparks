/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_ALIAS_SEARCH_H
#define SOLVE_ALIAS_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveAliasSearch : public Solve {
 public:
  SolveAliasSearch(class SPK *, int, char **);
  ~SolveAliasSearch();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  int allocated;
  int nevents;

  double *prob;
  double *p;
  double *q;
  int *j;
  int *hilo;
  int sk, sl, k, l;

  double sum;

  void free_arrays();
  void build_alias_table(int, double *);
  void table_dump(int);
  void check_table_consistency();

  int sort_hilo();
};

}

#endif
