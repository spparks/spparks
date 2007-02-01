/* ----------------------------------------------------------------------

The method sets up data structures that permit variate sampling from a
dynamic distribution with frequent updates. While the single variate 
generation time is logarithmic in distribution range here, 
the update time is constant and independent of distribution size.  
The distribution has to be bounded.
------------------------------------------------------------------------- */

#include "math.h"
#include <iostream>
#include "groups.h"
#include "random_park.h"

/* ---------------------------------------------------------------------- */
using namespace std;
using namespace SPPARKS;

Groups::Groups(double lo_in, double hi_in, int seed_in)
{
  my_group = NULL;
  my_group_i = NULL;
  group = NULL;
  group_size = NULL;
  i_group = NULL;
  group_hi = NULL;
  group_sum = NULL;
  empty_groups = NULL;

  random = new RandomPark(seed_in);
  hi = hi_in;
  lo = lo_in;
}
/* ---------------------------------------------------------------------- */
Groups::~Groups()
{
  release_group_space();
}
/* ----------------------------------------------------------------------
   Set up the partition groups data structures
   ------------------------------------------------------------------------- */
void Groups::partition_init(double *p, int size_in, int max_size_in)
{
  int i;
  range = hi - lo;

  size = size_in;
  max_size = max_size_in;

  overlg2 = 1.0/log(2.0); // set the value of the constant
  ngroups = 0;
  psum = 0.0;

  release_group_space();
  //set number of groups

  ngroups = static_cast<int> (log((double)max_size)*overlg2) + 1;
  allocate_group_space(ngroups);

  //init update data structures to produce error if not set
  for(i = 0; i < size; i++){
    my_group[i] = -1;
    my_group_i[i] = -1;
  }
  //  cout << "Created "<< ngroups << " groups."<<endl;

  int m, g;
  //calculate and allocate initial group storage

  //set group upper bounds and initialize empty/full vector
  double frac = range;
  for (g = 0; g < ngroups; g++){
    group_size[g] = 0;
    group_hi[g] = lo + frac;
    frac /= 2.0;
    empty_groups[g]=0;
  }
  //calculate initial group membership and initial sum
  for(i=0;i<size;i++){
    int gr = -static_cast<int>(log((p[i]-lo)/range)*overlg2);
    if (gr > ngroups-1) {
      gr = ngroups - 1;
      //      cout << "Distribution not flat in group "<<gr<<"."<<endl;
    }
    group_size[gr] ++;
    empty_groups[gr]=1;
    psum += p[i];
  }

  //set initial group size and allocate group storage
  for (g = 0; g < ngroups; g++){
    m = static_cast<int>(static_cast<double>(size)/pow(2.0,g));
    if(group_size[g]*2>m) m = group_size[g]*2;
    group[g] = new int[m];  
    group_size[g] = m;
  }
  //initialize group storage
  for (g = 0; g < ngroups; g++){
    i_group[g] = 0;
    group_sum[g] = 0.0;
    for(m = 0; m < group_size[g]; m++)
      group[g][m]=0;
  }
  //partition the initial distribution
  partition(p, lo, hi);
}
/* ----------------------------------------------------------------------
   Split distribution into a logarithmic number of groups by size
   ------------------------------------------------------------------------- */
void Groups::partition(double *p, double lo, double hi)
{
  double range = hi - lo;

  for(int j=0;j<size;j++){
    int g = -static_cast<int>(log((p[j]-lo)/range)*overlg2);
    if (g > ngroups-1) g = ngroups - 1; 
    group[g][i_group[g]] = j;
    my_group[j] = g;
    my_group_i[j] = i_group[g];
    i_group[g]++;
    group_sum[g] += p[j];
  }
}
/* ----------------------------------------------------------------------
   Change value of an element
   ------------------------------------------------------------------------- */
void Groups::alter_element(int j, double *p, double p_new)
{
  double p_old = p[j];
  double diff = p_new - p_old;

  //find new group membership
  int new_group = -static_cast<int>(log((p_new-lo)/range)*overlg2);
  if (new_group > ngroups-1) new_group = ngroups - 1; 

  //check if changed
  if(my_group[j] == new_group){
    //update group sum
    group_sum[new_group] += diff;
  }
  else{
    //remove from old group
    //put last element in place of removed
    int old_group = my_group[j];
    int old_group_i = my_group_i[j];
    int last = group[old_group][i_group[old_group]-1];

    group[old_group][old_group_i] = last;
    my_group_i[last] = old_group_i;
    i_group[old_group]--;
    group_sum[old_group] -= p_old;

    //add to new group
    if(i_group[new_group] >= group_size[new_group]-1) resize_group(new_group);
    group[new_group][i_group[new_group]] = j;
    my_group[j] = new_group;
    my_group_i[j] = i_group[new_group];
    i_group[new_group]++;
    group_sum[new_group] += p_new;
  }
  //update total sum
  psum += diff;
}
/* ----------------------------------------------------------------------
   Choose a group for new propensity element and add to group
   ------------------------------------------------------------------------- */
void Groups::add_element(int j, double *p)
{
  double p_in = p[j];
  //find group membership
  int g = -static_cast<int>(log((p_in-lo)/range)*overlg2);
  if (g > ngroups-1) g = ngroups - 1; 
  //  cout << "Adding propensity "<<p_in <<" to group "<<g<<endl;
  //check capacity and allocate space if needed
  if(i_group[g] >= group_size[g]-1) resize_group(g);
  if(size >= max_size) resize_inverse();
  //add to group
  group[g][i_group[g]] = j;
  my_group[j] = g;
  my_group_i[j] = i_group[g];
  i_group[g]++;
  group_sum[g] += p_in;
  
  //update total sum
  psum += p_in;
  //update size
  size ++;
}
/* ----------------------------------------------------------------------
   Extend group space
   ------------------------------------------------------------------------- */
void Groups::resize_group(int g)
{
  //double storage for group
  if(group_size[g] == 0) group_size[g] = 1;
  group_size[g] *= 2;
  int *tmpg = new int[group_size[g]];

  //copy contents
  for(int i=0;i<i_group[g]; i++) tmpg[i] = group[g][i];

  delete [] group[g];
  group[g] = tmpg;

   cout << "Resized group "<<g <<"  from "<<group_size[g]/2 <<" to "
       << group_size[g]<<"."<<endl;
}
/* ----------------------------------------------------------------------
   Extend inverse arrays
   ------------------------------------------------------------------------- */
void Groups::resize_inverse()
{
  //double storage for inverse
  if(max_size == 0) max_size = 1;
  max_size *= 2;
  int *tmp = new int[max_size];
  int *tmpi = new int[max_size];

  //copy contents
  for(int i=0;i<max_size; i++){
    tmp[i] = my_group[i];
    tmpi[i] = my_group_i[i];
  }

  delete [] my_group;
  delete [] my_group_i;

  my_group = tmp;
  my_group_i = tmpi;

//   cout << "Resized inverses from "<<max_size/2 <<" to "
//        << max_size<<"."<<endl;
}
/* ----------------------------------------------------------------------
   Sample distribution
   ------------------------------------------------------------------------- */
int Groups::sample(double *p)
{
  int r= -1;
  while(r<0){
    int grp = linear_select_group();
    //    cout << "group selected "<<grp <<endl;
    if(grp>-1) r =  sample_with_rejection(grp, p);
    //    cout << "group sampled "<<endl;
  }
  return r;
}
/* ----------------------------------------------------------------------
   Sample in-group distribution with rejection
   ------------------------------------------------------------------------- */
int Groups::sample_with_rejection(int g, double *p)
{
  int group_sample = static_cast<int>(i_group[g]*random->uniform());
  int sample = group[g][group_sample];
  
  if(p[sample] > group_hi[g]*random->uniform()) return sample; // accepted

  return -1;  //rejected
}
/* ----------------------------------------------------------------------
   Sample group sum distribution with partial sums
   ------------------------------------------------------------------------- */
int Groups::linear_select_group()
{
  int g;

  double compare = psum*random->uniform();
  double partial = 0.0;
  g = 0;

  while (1==1){
    partial += group_sum[g];
    if (partial > compare) return g;
    g++;
  }
  //  cout << "no selection of group"<<endl;
  return -1; //no group selected
}
/* ----------------------------------------------------------------------
   Allocate group arrays
   ------------------------------------------------------------------------- */
void Groups::allocate_group_space(int s)
{
  //allocate update data structures
  my_group = new int[max_size];
  my_group_i = new int[max_size];

  //allocate group data structures
  group_size = new int [s];
  i_group = new int [s];
  group_sum = new double [s];
  group_hi = new double [s];
  group = new int* [s];
  empty_groups = new int[s];
}
/* ----------------------------------------------------------------------
   Clear group arrays
   ------------------------------------------------------------------------- */
void Groups::release_group_space()
{
  if(my_group != NULL) delete [] my_group;
  if(my_group_i != NULL) delete [] my_group_i;
  if(group_size != NULL) delete [] group_size;
  if(i_group != NULL) delete [] i_group;
  if(group_sum != NULL) delete [] group_sum;
  if(group_hi != NULL) delete [] group_hi;
  if(empty_groups != NULL) delete [] empty_groups;
}

/* ----------------------------------------------------------------------
   Test sampling
   ------------------------------------------------------------------------- */
void Groups::test_sampling(double *p, int nsamples)
{
  int rjct = 0;

  int *count = new int[max_size];

  for(int i = 0; i < max_size; i++) count[i] = 0;

  for (int s = 0; s < nsamples; s++){
    int grp = linear_select_group();
    //    cout << "group selected "<<grp <<endl;
    int r =  sample_with_rejection(grp, p);
    //    cout << "group sampled "<<endl;
    if(r>-1) count[r]++;
    else rjct ++;
    //    cout <<"Sample group : "<< grp<<endl;
    //    cout << "Sample : "<<r<<" group "<<my_group[r]<<" index "
    //	 <<my_group_i[r]<<endl;
  }

  double ndsamples = static_cast<double> (nsamples-rjct);
  double current_max = 0.0;
  int current_max_i;
  double sumd = 0.0;
  double sumsqd = 0.0;
  for(int j = 0; j< size; j++){
    double ratio = static_cast<double>(count[j])*psum/ndsamples;
    double comp = abs(100.0*(p[j] - ratio)/(ratio + p[j]));
    if(current_max < comp){
      current_max = comp;
      current_max_i = j;
    }
    sumd += p[j] - ratio;
    sumsqd += (p[j] - ratio) * (p[j] - ratio);
//     cout << j << "  "<< p[j] << "  " 
//     	 << ratio<< "  "
//     	 << 100.0*(p[j] - ratio)/(ratio + p[j])<< "%"
//     	 <<endl;
    
  }
  sumsqd /= static_cast<double>(size);
  sumd /= static_cast<double>(size);
  cout << "***********************************************************"<<endl;
  cout << "Final distribution size = "<<size<<"."<<endl;
  cout << "Distribution sum = "<<psum<<"."<<endl;
  cout<< "Largest deviation = "<< current_max <<"%."<<endl;
  cout <<"The value of propensity for the largest deviation = "
       <<p[current_max_i]<<"."<<endl;
  cout << "Standard deviation = "<< sqrt(sumsqd - sumd*sumd) <<"."<<endl;
  cout << "Rejected "<< 100.0*(double)rjct/(double)nsamples 
       << "% uniform random numbers."<<endl;
  cout << "***********************************************************"<<endl;
  delete [] count;
}
/* ----------------------------------------------------------------------
   Output groups to screen diagnostic
   ------------------------------------------------------------------------- */
void Groups::group_diagnostic(double *p)
{
  //diagnostic
  cout <<"%%%%%%%%%%%%%%% Group diagnostic dump %%%%%%%%%%%%%%%%%%%%%"<<endl;
  cout << "Number of groups = "<<ngroups<<"."<<endl;
  for (int g = 0; g < ngroups; g++){
    cout << endl <<"Group "<< g<< " current/max "
 	 <<i_group[g]<<" / "<<group_size[g]<<" sums to " 
 	 <<group_sum[g]<<" with upper bound of "<<group_hi[g]
 	 <<"."<<endl<<endl;

    for (int m = 0; m < i_group[g]; m++){
      int l = group[g][m];
      cout << p[group[g][m]]<<"  "<<my_group[l]<<"  "
	   <<my_group_i[l]<<endl;

    }
    cout << endl 
	 <<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
  } 
  for(int i = 0; i< size; i++){
    //    cout << i <<"  " <<p[i]<<"  " <<p[i] - lo<< "  "<< my_group[i] <<"  "
    //     	 << my_group_i[i]<<endl;
    if(group[my_group[i]][my_group_i[i]] != i){
      cout << "At "<<i<<" the inverse table failed!"<<endl;
      cout <<" Producing "<< group[my_group[i]][my_group_i[i]]<<"."<<endl;
    }
  }
  
  cout << "Done group diagnostic .."<<endl;
}
/* ----------------------------------------------------------------------
   ------------------------------------------------------------------------- */
