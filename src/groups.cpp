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

  g2p = NULL;
  gcount = NULL;
  gmaxsize = NULL;
  ghibound = NULL;
  gpsum = NULL;

  p2g = NULL;
  p2g_index = NULL;
}

/* ---------------------------------------------------------------------- */

Groups::~Groups()
{
  release_memory();
  delete random;
}

/* ----------------------------------------------------------------------
   setup the group data structures
   each input propensity is between lo,hi inclusive or is zero
------------------------------------------------------------------------- */

void Groups::partition(double *p, int size_in)
{
  int g;

  size = size_in;
  psum = 0.0;

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
      ghibound[g] = lo + (g+1)*binsize;
    ghibound[ngroups-1] = hi;
  } else {
    double top = hi;
    for (int g = 0; g < ngroups; g++) {
      ghibound[g] = top;
      top *= 0.5;
    }
  }

  // count # of propensities in each group
  // psum = total propensity

  for (g = 0; g <= ngroups; g++) gcount[g] = 0;

  for (int i = 0; i < size; i++) {
    if (p[i] == 0.0) g = ngroups;
    else if (ngroups_flag) g = static_cast<int> ((p[i]-lo)*invbinsize);
    else {
      frexp(p[i]/hi,&g);
      g = -g;
      if (g < 0) g = 0;
    }

    gcount[g]++;
    psum += p[i];
  }
  
  // allocate per-group lists
  // ave = average group size times 1.5 factor

  int ave = static_cast<int> (1.5*size/ngroups);
  if (ngroups == 1) ave = size;
  if (ave == 0) ave = 1;

  for (g = 0; g <= ngroups; g++) {
    gmaxsize[g] = ave;
    if (gcount[g] > ave) gmaxsize[g] = gcount[g];
    memory->create(g2p[g],gmaxsize[g],"group:g2p");
  }

  // add index of each propensity to its group
  // each propensity is guaranteed to be 0.0 or between lo,hi inclusive

  for (g = 0; g <= ngroups; g++) {
    gcount[g] = 0;
    gpsum[g] = 0.0;
  }
  
  for (int i = 0; i < size; i++) {
    if (p[i] == 0.0) g = ngroups;
    else if (ngroups_flag) g = static_cast<int> ((p[i]-lo)*invbinsize);
    else {
      frexp(p[i]/hi,&g);
      g = -g;
      if (g < 0) g = 0;
    }

    p2g[i] = g;
    p2g_index[i] = gcount[g];
    g2p[g][gcount[g]++] = i;
    gpsum[g] += p[i];
  }

  // debug check

  //sanity_check(p);
}

/* ----------------------------------------------------------------------
   change value of a single element
   new propensity is between lo,hi inclusive or is zero
------------------------------------------------------------------------- */

void Groups::alter_element(int n, double *p, double p_new)
{
  double p_old = p[n];
  double diff = p_new - p_old;
  psum += diff;

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

  if (p2g[n] == new_group) gpsum[new_group] += diff;
  else {
    int old_group = p2g[n];
    int index = p2g_index[n];
    int lastprop = g2p[old_group][gcount[old_group]-1];

    g2p[old_group][index] = lastprop;
    p2g_index[lastprop] = index;
    gcount[old_group]--;
    gpsum[old_group] -= p_old;
    
    if (gcount[new_group] == gmaxsize[new_group])
      grow_group(new_group);

    p2g[n] = new_group;
    p2g_index[n] = gcount[new_group];
    g2p[new_group][gcount[new_group]++] = n;
    gpsum[new_group] += p_new;
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
  double thresh = psum*random->uniform();
  double partial = 0.0;

  int g;
  for (g = 0; g < ngroups; g++) {
    partial += gpsum[g];
    if (partial > thresh) return g;
  }
  return ngroups-1;
}

/* ----------------------------------------------------------------------
   sample in-group distribution with rejection
------------------------------------------------------------------------- */

int Groups::sample_with_rejection(int g, double *p)
{
  int index = static_cast<int>(gcount[g]*random->uniform());
  int sample = g2p[g][index];
  if (p[sample] > ghibound[g]*random->uniform()) return sample;
  return -1;
}

/* ----------------------------------------------------------------------
   double a group's size
------------------------------------------------------------------------- */

void Groups::grow_group(int g)
{
  gmaxsize[g] *= 2;
  memory->grow(g2p[g],gmaxsize[g],"group:g2p");
}

/* ----------------------------------------------------------------------
   allocate arrays for N+1 groups and for size propensities
   N = # of groups with non-zero propensities
   allocate one extra to store zero propensities
------------------------------------------------------------------------- */

void Groups::allocate_memory(int n)
{
  g2p = new int*[n+1];

  gcount = new int[n+1];
  gmaxsize = new int[n+1];
  gpsum = new double[n+1];
  ghibound = new double[n+1];

  memory->create(p2g,size,"group:p2g");
  memory->create(p2g_index,size,"group:p2g_index");
}

/* ----------------------------------------------------------------------
   deallocate group and propensity arrays
------------------------------------------------------------------------- */

void Groups::release_memory()
{
  if (g2p) {
    for (int g = 0; g <= ngroups; g++) memory->destroy(g2p[g]);
    delete [] g2p;
    g2p = NULL;
  }

  delete [] gcount;
  delete [] gmaxsize;
  delete [] gpsum;
  delete [] ghibound;

  memory->destroy(p2g);
  memory->destroy(p2g_index);
}

/* ----------------------------------------------------------------------
   sanity check that group data structure is consistent
------------------------------------------------------------------------- */

void Groups::sanity_check(double *p)
{
  // check that per-group count and sum = totals within EPSILON

  int count = 0;
  double sum = 0.0;

  for (int g = 0; g <= ngroups; g++) {
    count += gcount[g];
    sum += gpsum[g];
  }

  if (count != size) printf("Bad group count: %d %d\n",count,size);
  if (fabs(sum-psum) > EPSILON) printf("Bad group sum: %g %g\n",sum,psum);

  // check that all propensities in group are within group bounds

  int flag = 0;
  for (int g = 0; g <= ngroups; g++) {
    double glo,ghi;
    if (g < ngroups) {
      if (ngroups_flag) {
	if (g == 0) glo = lo;
	else glo = ghibound[g-1];
	ghi = ghibound[g];
      } else {
	if (g < ngroups-1) glo = ghibound[g+1];
	else glo = lo;
	ghi = ghibound[g];
      }
    } else glo = ghi = 0.0;

    for (int i = 0; i < gcount[g]; i++) {
      double pone = p[g2p[g][i]];
      if (pone < glo || pone > ghi) {
	printf("Group mismatch: %d %d %d: %g %g %g\n",
	       g,i,g2p[g][i],glo,pone,ghi);
	flag++;
      }
    }
  }

  if (flag) printf("Mis-match of propensity with group: %d\n",flag);
}
