/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef SolveClass
SolveStyle(tree,SolveTree)

#else

  //#ifndef SPK_SOLVE_TREE_H
  //#define SPK_SOLVE_TREE_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveTree : public Solve {
 public:
  SolveTree(class SPPARKS *, int, char **);
  ~SolveTree();
  SolveTree *clone();

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
  class RandomPark *random;
  int nevents;
  int ntotal;

  double *tree;
  int offset;
  int allocated;
};

}

#endif
