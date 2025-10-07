/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_GROUPS_H
#define SPK_GROUPS_H

#include "pointers.h"

namespace SPPARKS_NS {

class Groups : protected Pointers {
 public:
  Groups(class SPPARKS *, double, double, int);
  ~Groups();

  void partition(double *,int);
  void alter_element(int, double *, double);
  int sample(double *);

 private:
  int size;             // number of propensities
  double psum;          // sum of all propensities
  double hi,lo;         // inclusive propensity bounds
  double invbinsize;    // inverse width of propensity bin when user-specified

  int ngroups;          // number of groups
  int ngroups_flag;     // 0 for logarithmic, N for N equal-spaced groups
  int **g2p;            // list of propensity indices in each group
  int *gcount;          // # of propensities in each group
  int *gmaxsize;        // max # of propensities each group is allocated for
  double *ghibound;     // upper bound of propensity for each group
  double *gpsum;        // sum of propensities for each group

  int *p2g;             // p2g[n] = which group propensity N is assigned to
  int *p2g_index;       // p2gindex[n] = index within group of propensity N

  class RandomPark *random;

  int linear_select_group();
  int sample_with_rejection(int, double *);
  void grow_group(int);
  void allocate_memory(int);
  void release_memory();
  int find_group(double pone);
  void sanity_check(double *p);
};

}

#endif
