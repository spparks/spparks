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

#include "math.h"
#include "stdlib.h"
#include "groups.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

Groups::Groups(SPPARKS *spk, double hi_in, double lo_in, int ng_in) : 
  Pointers(spk)
{
  random = new RandomPark(ranmaster->uniform());

  hi = hi_in;
  lo = lo_in;
  ngroups = ng_in;
  if (ngroups == 0) ngroups_flag = 0;
  else ngroups_flag = 1;

  group = NULL;
  group_maxsize = NULL;
  group_size = NULL;
  group_hi = NULL;
  group_sum = NULL;

  my_group = NULL;
  my_group_i = NULL;
}

/* ---------------------------------------------------------------------- */

Groups::~Groups()
{
  release_memory();
  delete random;
}

/* ----------------------------------------------------------------------
   setup the group data structures
------------------------------------------------------------------------- */

void Groups::partition(double *p, int size_in)
{
  int g;

  size = size_in;
  sum = 0.0;

  // for ngroups_flag = 1:
  //   propensity = lo will be in group 0
  //   propensity = hi will be in ngroup-1
  // for ngroups_flag = 0:
  //   frexp(hi/hi) will be in group 0
  //   frexp(lo/hi) will be in ngroup-1
  // in both cases propensity 0.0 will be put in ngroup
  // allocate_group_space thus allocates ngroup+1 length arrays 

  release_memory();
  if (ngroups_flag == 0) {
    frexp(lo/hi,&ngroups);
    ngroups = -ngroups;
    ngroups++;
  }
  allocate_memory(ngroups);

  // compute upper propensity bound for each group
  // except group that stores zero propensities
  // add epsilon to binsize to insure propensity = hi falls in last group

  if (ngroups_flag) {
    double binsize = (hi-lo)/(ngroups-EPSILON);
    invbinsize = 1.0/binsize;
    for (int g = 0; g < ngroups; g++)
      group_hi[g] = lo + (g+1)*binsize;
    group_hi[ngroups-1] = hi;
  } else {
    double top = hi;
    for (int g = 0; g < ngroups; g++) {
      group_hi[g] = top;
      top *= 0.5;
    }
  }

  // count # of propensities in each group
  // sum = total propensity

  for (g = 0; g <= ngroups; g++) group_size[g] = 0;

  for (int i = 0; i < size; i++) {
    if (p[i] == 0.0) g = ngroups;
    else if (ngroups_flag) g = static_cast<int> ((p[i]-lo)*invbinsize);
    else {
      frexp(p[i]/hi,&g);
      g = -g;
      if (g < 0) g = 0;
    }

    group_size[g]++;
    sum += p[i];
  }
  
  // allocate per-group lists
  // ave = average group size times 1.5 factor

  int ave = static_cast<int> (1.5*size/ngroups);
  if (ngroups == 1) ave = size;
  if (ave == 0) ave = 1;

  for (g = 0; g <= ngroups; g++) {
    group_maxsize[g] = ave;
    if (group_size[g] > ave) group_maxsize[g] = group_size[g];
    group[g] = (int *) memory->smalloc(group_maxsize[g]*sizeof(int),
				       "group:group");
  }

  // add index of each propensity to its group
  // each propensity is guaranteed to be 0.0 or between lo,hi inclusive

  for (g = 0; g <= ngroups; g++) {
    group_size[g] = 0;
    group_sum[g] = 0.0;
  }
  
  for (int i = 0; i < size; i++) {
    if (p[i] == 0.0) g = ngroups;
    else if (ngroups_flag) g = static_cast<int> ((p[i]-lo)*invbinsize);
    else {
      frexp(p[i]/hi,&g);
      g = -g;
      if (g < 0) g = 0;
    }

    my_group[i] = g;
    my_group_i[i] = group_size[g];
    group[g][group_size[g]++] = i;
    group_sum[g] += p[i];
  }

  // debug check

  //sanity_check(p);
}

/* ----------------------------------------------------------------------
   change value of a single element
------------------------------------------------------------------------- */

void Groups::alter_element(int n, double *p, double p_new)
{
  double p_old = p[n];
  double diff = p_new - p_old;
  sum += diff;

  // compute new group for p_new

  int new_group;
  if (p_new == 0.0) new_group = ngroups;
  else if (ngroups_flag) new_group = static_cast<int> ((p_new-lo)*invbinsize);
  else {
    frexp(p_new/hi,&new_group);
    new_group = -new_group;
    if (new_group < 0) new_group = 0;
  }
  
  // if group changed, delete propensity from old group, add to new group
  // index = where in old group this propensity is
  // lastprop = which propensity is last in old group

  if (my_group[n] == new_group) group_sum[new_group] += diff;
  else {
    int old_group = my_group[n];
    int index = my_group_i[n];
    int lastprop = group[old_group][group_size[old_group]-1];

    group[old_group][index] = lastprop;
    my_group_i[lastprop] = index;
    group_size[old_group]--;
    group_sum[old_group] -= p_old;
    
    if (group_size[new_group] == group_maxsize[new_group])
      grow_group(new_group);

    my_group[n] = new_group;
    my_group_i[n] = group_size[new_group];
    group[new_group][group_size[new_group]++] = n;
    group_sum[new_group] += p_new;
  }

  // debug check

  //p[n] = p_new;
  //sanity_check(p);
}

/* ----------------------------------------------------------------------
   sample distribution
------------------------------------------------------------------------- */

int Groups::sample(double *p)
{
  int g;

  int r = -1;
  while (r < 0) {
    g = linear_select_group();
    r = sample_with_rejection(g,p);
  }
  return r;
}

/* ----------------------------------------------------------------------
   sample group sum distribution in linear fashion
------------------------------------------------------------------------- */

int Groups::linear_select_group()
{
  double thresh = sum*random->uniform();
  double partial = 0.0;

  int g;
  for (g = 0; g < ngroups; g++) {
    partial += group_sum[g];
    if (partial > thresh) return g;
  }
  return ngroups-1;
}

/* ----------------------------------------------------------------------
   sample in-group distribution with rejection
------------------------------------------------------------------------- */

int Groups::sample_with_rejection(int g, double *p)
{
  int index = static_cast<int>(group_size[g]*random->uniform());
  int sample = group[g][index];
  if (p[sample] > group_hi[g]*random->uniform()) return sample;
  return -1;
}

/* ----------------------------------------------------------------------
   double a group's size
------------------------------------------------------------------------- */

void Groups::grow_group(int g)
{
  group_maxsize[g] *= 2;
  group[g] = (int *) memory->srealloc(group[g],group_maxsize[g]*sizeof(int),
				      "group:group");
}

/* ----------------------------------------------------------------------
   allocate arrays for N+1 groups and for size propensities
   N = # of groups with non-zero propensities
   allocate one extra to store zero propensities
------------------------------------------------------------------------- */

void Groups::allocate_memory(int n)
{
  group = new int*[n+1];

  group_maxsize = new int [n+1];
  group_size = new int [n+1];
  group_sum = new double [n+1];
  group_hi = new double [n+1];

  my_group = (int *) memory->smalloc(size*sizeof(int),"group:my_group");
  my_group_i = (int *) memory->smalloc(size*sizeof(int),"group:my_group_i");
}

/* ----------------------------------------------------------------------
   deallocate group and propensity arrays
------------------------------------------------------------------------- */

void Groups::release_memory()
{
  if (group) {
    for (int g = 0; g <= ngroups; g++) memory->sfree(group[g]);
    delete [] group;
    group = NULL;
  }

  delete [] group_maxsize;
  delete [] group_size;
  delete [] group_sum;
  delete [] group_hi;

  memory->sfree(my_group);
  memory->sfree(my_group_i);
}

/* ----------------------------------------------------------------------
   sanity check that group data structure is consistent
------------------------------------------------------------------------- */

void Groups::sanity_check(double *p)
{
  // check that per-group count and sum = totals within EPSILON

  int gcount = 0;
  double gsum = 0.0;

  for (int g = 0; g <= ngroups; g++) {
    gcount += group_size[g];
    gsum += group_sum[g];
  }

  if (gcount != size) printf("Bad group count: %d %d\n",gcount,size);
  if (fabs(gsum-sum) > EPSILON) printf("Bad group sum: %g %g\n",gsum,sum);

  // check that all propensities in group are within group bounds

  int flag = 0;
  for (int g = 0; g <= ngroups; g++) {
    double glo,ghi;
    if (g < ngroups) {
      if (ngroups_flag) {
	if (g == 0) glo = lo;
	else glo = group_hi[g-1];
	ghi = group_hi[g];
      } else {
	if (g < ngroups-1) glo = group_hi[g+1];
	else glo = lo;
	ghi = group_hi[g];
      }
    } else glo = ghi = 0.0;

    for (int i = 0; i < group_size[g]; i++) {
      double pone = p[group[g][i]];
      if (pone < glo || pone > ghi) {
	printf("Group mismatch: %d %d %d: %g %g %g\n",
	       g,i,group[g][i],glo,pone,ghi);
	flag++;
      }
    }
  }

  if (flag) printf("Mis-match of propensity with group: %d\n",flag);
}
