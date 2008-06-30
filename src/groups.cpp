/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "groups.h"
#include "random_park.h"
#include "memory.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Groups::Groups(SPPARKS *spk, double lo_in, double hi_in, int seed_in, 
		 int ng_flag, int ngr_in) : SysPtr(spk)
{
  my_group = NULL;
  my_group_i = NULL;
  group = NULL;
  group_size = NULL;
  i_group = NULL;
  group_hi = NULL;
  group_sum = NULL;
  empty_groups = NULL;
  nempty = 0;

  random = new RandomPark(seed_in);
  hi = hi_in;
  lo = lo_in;
  ngroups_flag = ng_flag;
  ngroups = ngr_in;
}

/* ---------------------------------------------------------------------- */

Groups::~Groups()
{
  for (int g = 0; g < ngroups; g++) memory->sfree(group[g]);
  delete [] group;
  release_group_space();
  delete random;
}

/* ----------------------------------------------------------------------
   set up the groups data structures
------------------------------------------------------------------------- */

void Groups::partition_init(double *p, int size_in)
{
  int i;
  range = hi/lo;

  max_size = size = size_in;

  overlg2 = 1.0/log(2.0);

  if (!ngroups_flag) ngroups = 0;
  psum = 0.0;

  release_group_space();
  if (!ngroups_flag) ngroups = static_cast<int> (log(range)*overlg2) + 1;
  allocate_group_space(ngroups);

  // init update data structures to produce error if not set

  for (i = 0; i < size; i++) {
    my_group[i] = -1;
    my_group_i[i] = -1;
  }

  int m,g;

  // calculate and allocate initial group storage
  // set group upper bounds and initialize empty/full vector

  double frac;
  if (ngroups_flag)
    frac = hi/static_cast<double>(ngroups);
  else
    frac = hi;

  for (g = 0; g < ngroups; g++){
    group_size[g] = 0;
    empty_groups[g]=0;
    if (ngroups_flag)
      group_hi[g] = (g+1)*frac;
    else{
      group_hi[g] = frac;
      frac /= 2.0;
    }
  }

  // calculate initial group membership and initial sum

  int gr;
  for (i = 0; i < size; i++) {
    if (ngroups_flag) gr = static_cast<int>(p[i]/frac);
    else if (p[i] > 1.0e-20) gr = -static_cast<int>(log(p[i]/hi)*overlg2);
    else gr = ngroups - 1;
    if (gr > ngroups-1) gr = ngroups - 1;

    group_size[gr] ++;
    empty_groups[gr]=0;
    psum += p[i];
  }
  
  // set initial group size and allocate group storage

  for (g = 0; g < ngroups; g++) {
    if (ngroups_flag) m = static_cast<int>
      (1.5*static_cast<double>(size)/static_cast<double>(ngroups));
    else m = static_cast<int>(static_cast<double>(size)/pow(2.0,g));
    if (group_size[g] > m) m = group_size[g];
    group[g] = (int *) memory->smalloc(m*sizeof(int),"group:group");
    group_size[g] = m;
  }

  // initialize group storage

  for (g = 0; g < ngroups; g++){
    i_group[g] = 0;
    group_sum[g] = 0.0;
    for (m = 0; m < group_size[g]; m++) group[g][m] = 0;
  }

  // partition the initial distribution

  partition(p, lo, hi);
}

/* ----------------------------------------------------------------------
   split distribution into a number of groups by size
------------------------------------------------------------------------- */

void Groups::partition(double *p, double lo, double hi)
{
  double range = hi - lo;
  int g = 0;
  double tst;

  // equal fragments

  if (ngroups_flag) {
    frac = hi/static_cast<double>(ngroups);
    for(int j=0;j<size;j++){
      g = static_cast<int>(p[j]/frac);
      group[g][i_group[g]] = j;
      my_group[j] = g;
      my_group_i[j] = i_group[g];
      i_group[g]++;
      group_sum[g] += p[j];
    }

  // logarithmic fragments

  } else {
    frac = hi;
    for(int j=0;j<size;j++){
      if (p[j] > 1.0e-20) {
	tst = frexp(p[j]/hi,&g);
	g = -g;
      } else g = ngroups - 1;
      if (g > ngroups-1) g = ngroups - 1; 
      group[g][i_group[g]] = j;
      my_group[j] = g;
      my_group_i[j] = i_group[g];
      i_group[g]++;
      group_sum[g] += p[j];
      
    }
  }

  for (int g = 0; g < ngroups; g++) 
    if(i_group[g]>0) {
      empty_groups[g] = 1;
      nempty ++;
    }
}

/* ----------------------------------------------------------------------
   change value of a single element
------------------------------------------------------------------------- */

void Groups::alter_element(int j, double *p, double p_new)
{
  double p_old = p[j];
  double diff = p_new - p_old;
  psum += diff;

  // find new group membership

  int new_group,old_group;

  if (ngroups_flag) new_group = static_cast<int>(p_new/frac);
  else {
    if (p_new > lo) {
      double tmp = frexp(p[j]/hi,&new_group); 
      new_group = -new_group;
    } else new_group = ngroups - 1; 
  }
  
  // check if changed

  if (my_group[j] == new_group) group_sum[new_group] += diff;
  else {
    old_group = my_group[j];
    int old_group_i = my_group_i[j];
    int last = group[old_group][i_group[old_group]-1];
    
    group[old_group][old_group_i] = last;
    my_group_i[last] = old_group_i;
    i_group[old_group]--;
    group_sum[old_group] -= p_old;
    
    if (i_group[new_group] > group_size[new_group]-1) grow_group(new_group);

    group[new_group][i_group[new_group]] = j;
    my_group[j] = new_group;
    my_group_i[j] = i_group[new_group];
    i_group[new_group]++;
    group_sum[new_group] += p_new;
  }
}

/* ----------------------------------------------------------------------
   sample distribution
------------------------------------------------------------------------- */

int Groups::sample(double *p)
{
  int grp;
  int cnt = 0;

  int r = -1;
  while (r < 0) {
    grp = linear_select_group();
    if (grp > -1) r = sample_with_rejection(grp,p);
  }
  return r;
}

/* ----------------------------------------------------------------------
   double a group's size
------------------------------------------------------------------------- */

void Groups::grow_group(int g)
{
  if (group_size[g] == 0) group_size[g] = 1;
  group_size[g] *= 2;
  group[g] = (int *) memory->srealloc(group[g],group_size[g]*sizeof(int),
				      "group:group");
}

/* ----------------------------------------------------------------------
   sample in-group distribution with rejection
------------------------------------------------------------------------- */

int Groups::sample_with_rejection(int g, double *p)
{
  int group_sample = static_cast<int>(i_group[g]*random->uniform());
  int sample = group[g][group_sample];
  
  if (p[sample] > group_hi[g]*random->uniform()) return sample;
  return -1;
}

/* ----------------------------------------------------------------------
   sample group sum distribution with partial sums
------------------------------------------------------------------------- */

int Groups::linear_select_group()
{
  int g;

  double compare = psum*random->uniform();
  double partial = 0.0;
  g = 0;

  while (g < ngroups+1) {
    partial += group_sum[g];
    if (partial > compare) return g;
    g++;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   allocate group arrays
------------------------------------------------------------------------- */

void Groups::allocate_group_space(int s)
{
  my_group = new int[max_size];
  my_group_i = new int[max_size];

  group_size = new int [s];
  i_group = new int [s];
  group_sum = new double [s];
  group_hi = new double [s];
  group = new int*[s];
  empty_groups = new int[s];
}

/* ----------------------------------------------------------------------
   clear group arrays
------------------------------------------------------------------------- */

void Groups::release_group_space()
{
  if (my_group != NULL) delete [] my_group;
  if (my_group_i != NULL) delete [] my_group_i;
  if (group_size != NULL) delete [] group_size;
  if (i_group != NULL) delete [] i_group;
  if (group_sum != NULL) delete [] group_sum;
  if (group_hi != NULL) delete [] group_hi;
  if (empty_groups != NULL) delete [] empty_groups;
}
