/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef GROUPS_H
#define GROUPS_H

namespace SPPARKS {

class Groups {
 public:
  Groups(double, double, int, bool, int);
  ~Groups();

  void partition_init(double *,int, int);
  void partition(double *, double, double);

  int sample(double *);
  void remove_element(int, double *);
  void alter_element(int, double *, double);
  void add_element(int, double *);

  void group_diagnostic(double *);
  void test_sampling(double *, int);
  int diag_cnt;

 private:
  double overlg2; // 1.0/log(2.0)

  int *my_group; //inverse group index
  int *my_group_i; //inverse group location index 
  int size; //current distribution size
  int max_size; //current allocated distribution storage
  double psum; //current distribution sum
  double hi, lo, range; //distribution parameters

  int ngroups; //number of groups
  bool ngroups_flag;
  int *group_size; //size of groups
  int **group; //group members
  int *i_group; //current group max index
  double *group_hi; //current group sup
  double *group_sum; // current group sum
  int *empty_groups;  //0 for empty, 1 for non-empty
  int nempty;
  double frac;

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
