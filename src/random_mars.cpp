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
#include "random_mars.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- 
    Marsaglia RNG
------------------------------------------------------------------------ */

RanMars::RanMars(SPPARKS *spk) : Pointers(spk)
{
  initflag = 0;
  u = NULL;
}

/* ---------------------------------------------------------------------- */

RanMars::~RanMars()
{
  delete [] u;
}

/* ---------------------------------------------------------------------- */

void RanMars::init(int seed)
{
  int ij,kl,i,j,k,l,ii,jj,m;
  double s,t;

  initflag = 1;

  // assume input seed is positive value > 0
  // insure seed is from 1 to 900,000,000 inclusive
  
  while (seed > 900000000) seed -= 900000000;

  u = new double[97+1];

  ij = (seed-1)/30082;
  kl = (seed-1) - 30082*ij;
  i = (ij/177) % 177 + 2;
  j = ij %177 + 2;
  k = (kl/169) % 178 + 1;
  l = kl % 169;
  for (ii = 1; ii <= 97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj = 1; jj <= 24; jj++) {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s = s + t;
      t = 0.5*t;
    }
    u[ii] = s;
  }
  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;
  uniform();
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RanMars::uniform()
{
  if (!initflag) error->all(FLERR,"Seed command has not been used");

  double uni = u[i97] - u[j97];
  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97 == 0) i97 = 97;
  j97--;
  if (j97 == 0) j97 = 97;
  c -= cd;
  if (c < 0.0) c += cm;
  uni -= c;
  if (uni < 0.0) uni += 1.0;
  return uni;
}
