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

#include "string.h"
#include "stdlib.h"
#include "potential.h"
#include "error.h"

#include "style_pair.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Potential::Potential(SPPARKS *spk) : Pointers(spk)
{
  pair_style = NULL;
  pair = NULL;
}

/* ---------------------------------------------------------------------- */

Potential::~Potential()
{
  delete [] pair_style;
  delete pair;
}

/* ---------------------------------------------------------------------- */

void Potential::init()
{
  if (pair) pair->init();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script
------------------------------------------------------------------------- */

void Potential::create_pair(const char *style)
{
  delete [] pair_style;
  if (pair) delete pair;

  pair = new_pair(style);
  int n = strlen(style) + 1;
  pair_style = new char[n];
  strcpy(pair_style,style);
}

/* ----------------------------------------------------------------------
   generate a pair class
------------------------------------------------------------------------- */

Pair *Potential::new_pair(const char *style)
{
  if (strcmp(style,"none") == 0) return NULL;

#define PAIR_CLASS
#define PairStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(spk);
#include "style_pair.h"
#undef PAIR_CLASS

  else error->all(FLERR,"Invalid pair style");
  return NULL;
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to nmax,
     (3) i* = i to nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Potential::bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),1);
    nhi = MIN(atoi(str),nmax);
  } else if (strlen(str) == 1) {
    nlo = 1;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 1;
    nhi = MIN(atoi(ptr+1),nmax);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),1);
    nhi = nmax;
  } else {
    nlo = MAX(atoi(str),1);
    nhi = MIN(atoi(ptr+1),nmax);
  }
}
