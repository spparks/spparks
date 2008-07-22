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
#include "random_ranecu.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

RandomRanecu::RandomRanecu(int seed_init1, int seed_init2)
{
  seed1 = seed_init1;
  seed2 = seed_init2;
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RandomRanecu::uniform()
{
  if (seed1 > 2147483562) seed1 /= 2;
  if (seed2 > 2147483398) seed2 /= 2;

  int k = seed1/53668;
  seed1 = 40014 * (seed1 - k*53668) - k*12211;
  if (seed1 < 0) seed1 += 2147483563;

  k = seed2/52774;
  seed2 = 40692 * (seed2 - k*52774) - k*3791;
  if (seed2 < 0) seed2 += 2147483399;

  int iz = seed1 - seed2;
  if (iz < 1) iz = iz + 2147483562;

  return iz*4.656613e-10;
}

/* ----------------------------------------------------------------------
   integer RN between 1 and N inclusive
------------------------------------------------------------------------- */

int RandomRanecu::irandom(int n)
{
  int i = (int) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}
