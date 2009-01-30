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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_tree.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

SolveTree::SolveTree(SPPARKS *spk, int narg, char **arg) : 
  Solve(spk, narg, arg)
{
  if (narg != 1) error->all("Illegal solve command");

  random = new RandomPark(ranmaster->uniform());
  allocated = 0;
  tree = NULL;
}

/* ---------------------------------------------------------------------- */

SolveTree::~SolveTree()
{
  free_arrays();
  delete random;
}

/* ---------------------------------------------------------------------- */

SolveTree *SolveTree::clone()
{
  int narg = 1;
  char *arg[1];
  arg[0] = style;

  SolveTree *ptr = new SolveTree(spk,narg,arg);

  return ptr;
}

/* ---------------------------------------------------------------------- */

void SolveTree::init(int n, double *propensity)
{
  ntotal = 0;
  offset = 0;

  // memory allocation

  if (allocated) free_arrays();
  allocated = 1;

  nevents = n;
  sum = 0;
  num_active = 0;

  // m = value such that 2^m >= nevents

  int m = 0;
  long long int neat = 1;

  while (neat < nevents) {neat *=2; m++;}

  // create tree of length ntotal
  // use 64bit math, so we can handle ntotal = 2^31-1

  ntotal = 2LL*neat - 1;
  tree = new double[ntotal];
  offset = neat - 1;

  for (int i = 0; i < ntotal; i++) tree[i] = 0.0;
  for (int i = offset; i < offset + n; i++) 
    tree[i] = propensity[i-offset];
  sum_tree();

}

/* ---------------------------------------------------------------------- */

void SolveTree::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) set(indices[i],propensity[indices[i]]);
}

/* ---------------------------------------------------------------------- */

void SolveTree::update(int n, double *propensity)
{
  set(n,propensity[n]);
}

/* ---------------------------------------------------------------------- */

void SolveTree::resize(int new_size, double *propensity)
{
  init(new_size,propensity);
}

/* ---------------------------------------------------------------------- */

int SolveTree::event(double *pdt)
{
  int m;
  double r2;

  if (sum == 0.0) return -1;

  r2 = random->uniform();
  m = find(r2*sum);
  
  *pdt = -1.0/sum * log(random->uniform());

  return m;
}

/* ----------------------------------------------------------------------
   sum entire tree, all nodes are computed
------------------------------------------------------------------------- */

void SolveTree::sum_tree()
{
  int child1,child2;
  for (int parent = offset-1; parent >= 0; parent--) {
    child1 = 2*parent + 1;
    child2 = 2*parent + 2;
    tree[parent] = tree[child1] + tree[child2];
  }

  // update total propensity

  sum = tree[0];

  // update number of active events

  num_active = 0;
  for (int i = offset; i < ntotal; i++) 
    if (tree[i] > 0.0) num_active++;
}

/* ----------------------------------------------------------------------
   set propensity[i] to value
   recompute sum tree for all its ancestors
------------------------------------------------------------------------- */

void SolveTree::set(int i, double value)
{
  int parent,sibling;

  // update number of active events

  if (tree[offset+i] > 0.0) num_active--;
  if (value > 0.0) num_active++;

  tree[offset+i] = value;

  // I walks tree from leaf to root, summing children at each step
  // left child is odd index, right child is even index

  i += offset;
  while (i > 0) {
    if (i % 2) sibling = i + 1;
    else sibling = i - 1;
    parent = (i-1)/2;
    tree[parent] = tree[i] + tree[sibling];
    i = parent;
  }

  // update total propensity

  sum = tree[0];
}

/* ----------------------------------------------------------------------
   value = uniform RN from 0 to tree[0]
   return index (0 to M-1) of propensity bin it falls in
------------------------------------------------------------------------- */

int SolveTree::find(double value)
{
  int i,leftchild;

  // I walks tree from root to appropriate leaf
  // value is modified when right branch of tree is traversed

  i = 0;
  while (i < offset) {
    leftchild = 2*i + 1;
    if (value <= tree[leftchild]) i = leftchild;
    else {
      value -= tree[leftchild];
      i = leftchild + 1;
    }
  }
  return i - offset;
}

/* ----------------------------------------------------------------------
   free arrays used by solver
------------------------------------------------------------------------- */

void SolveTree::free_arrays()
{
  delete [] tree;
}
