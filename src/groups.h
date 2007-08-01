/* ----------------------------------------------------------------------

The method sets up data structures that permit variate sampling from a
dynamic distribution with frequent updates. While the single variate 
generation time is logarithmic in distribution range here, 
the update time is constant and independent of distribution size.  
The distribution has to be bounded.
------------------------------------------------------------------------- */

#ifndef GROUPS_H
#define GROUPS_H

#include "random_park.h"
#include <cmath>
#include <iostream>

using namespace std;
namespace SPPARKS {
class Groups {
 public:

  Groups(double, double, int, bool, int);
  ~Groups();

  //setup
  void partition_init(double *,int, int);
  void partition(double *, double, double);

  //updates
  void remove_element(int, double *); // remove a propensity 
  void alter_element(int, double *, double); //alter a propensity
  void add_element(int, double *); // add a propensity


  //diagnostics
  void group_diagnostic(double *); // output groups to screen
  void test_sampling(double *, int);   //test sampling statistics 
                                  //against distribution
  int sample(double *); //draw a variate from current distribution
 private:
  //constants
  double overlg2; // 1.0/log(2.0)
  //group parameters and storage

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

  void allocate_group_space(int);
  void release_group_space();
  void resize_group(int);
  void resize_inverse();

  //sampling
  class RandomPark *random;
  int sample_with_rejection(int, double *);
  int linear_select_group();


};

}

#endif

