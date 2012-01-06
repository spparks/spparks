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

#include "spktype.h"
#include "math.h"
#include "random_park.h"

using namespace SPPARKS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ---------------------------------------------------------------------- 
   Park/Miller RNG
   assume iseed is a positive int
------------------------------------------------------------------------ */

RandomPark::RandomPark(int iseed)
{
  seed = iseed;
}

/* ---------------------------------------------------------------------- 
   set seed to positive int
   assume 0.0 <= rseed < 1.0
------------------------------------------------------------------------ */

RandomPark::RandomPark(double rseed)
{
  seed = static_cast<int> (rseed*IM);
  if (seed == 0) seed = 1;
}

/* ---------------------------------------------------------------------- 
   reset seed to a positive int based on rseed and offset
   assume 0.0 <= rseed < 1.0 and offset is an int >= 0
   fmod() insures no overflow when static cast to int
   warmup the new RNG if requested
   typically used to setup one RN generator per proc or site or particle
------------------------------------------------------------------------ */

void RandomPark::reset(double rseed, int offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  for (int i = 0; i < warmup; i++) uniform();
}

void RandomPark::tagreset(double rseed, tagint offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  for (int i = 0; i < warmup; i++) uniform();
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RandomPark::uniform()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   integer RN between 1 and N inclusive
------------------------------------------------------------------------- */

int RandomPark::irandom(int n)
{
  int i = (int) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}

/* ----------------------------------------------------------------------
   tagint RN between 1 and N inclusive
------------------------------------------------------------------------- */

tagint RandomPark::tagrandom(tagint n)
{
  tagint i = (tagint) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}

/* ----------------------------------------------------------------------
   bigint RN between 1 and N inclusive
------------------------------------------------------------------------- */

bigint RandomPark::bigrandom(bigint n)
{
  bigint i = (bigint) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}
