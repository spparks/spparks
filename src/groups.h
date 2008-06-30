/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef GROUPS_H
#define GROUPS_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Groups : protected SysPtr {
 public:
  Groups(class SPPARKS *, double, double, int, int, int);
  ~Groups();

  void partition_init(double *,int);
  void alter_element(int, double *, double);
  int sample(double *);

 private:
  double overlg2;

  int *my_group;   // inverse group index
  int *my_group_i; // inverse group location index 
  int size;        // current distribution size
  int max_size;    // current allocated distribution storage
  double psum;     // current distribution sum
  double hi,lo,range;  // distribution parameters

  int ngroups;          // number of groups
  int ngroups_flag;
  int *group_size;      // size of groups
  int **group;          // group members
  int *i_group;         // current group max index
  double *group_hi;     // current group sup
  double *group_sum;    // current group sum
  int *empty_groups;    // 0 for empty, 1 for non-empty
  int nempty;
  double frac;

  void partition(double *, double, double);
  void remove_element(int, double *);
  void add_element(int, double *);

  void allocate_group_space(int);
  void release_group_space();
  void grow_group(int);
  void shrink_group(int);
  void resize_group(int);
  void resize_inverse();

  class RandomPark *random;
  int sample_with_rejection(int, double *);
  int linear_select_group();
};

}

#endif
