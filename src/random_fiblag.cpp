/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "random_fiblag.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

RandomFibLag::RandomFibLag(int seed_init)
{
  seed = seed_init;
  if (seed > 900000000) seed = (int) 0.41*seed;

  int ij = seed/30082;
  int kl = seed - 30082*ij;
  int i = (ij/177) % 177 + 2;
  int j = ij % 177 + 2;
  int k = (kl/169) % 178 + 1;
  int l = kl % 169;

  double s,t;
  int ii,jj,m;
  for (ii = 0; ii < 56; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj = 1; jj <= 24; jj++) {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s += t;
      t *= 0.5;
    }
    fib[ii] = s - (int) s;
  }

  i0 = 0;
  i8 = 8;
  i16 = 16;
  i24 = 24;
  i55 = 55;
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RandomFibLag::uniform()
{
  fib[i0] = fib[i8] + fib[i16] + fib[i24] + fib[i55];
  fib[i0] -= (int) fib[i0];
  double random = fib[i0];
  seed = static_cast<int> (2147483647*random);

  i0--;
  if (i0 < 0) i0 = 55;
  i8--;
  if (i8 < 0) i8 = 55;
  i16--;
  if (i16 < 0) i16 = 55;
  i24--;
  if (i24 < 0) i24 = 55;
  i55--;
  if (i55 < 0) i55 = 55;

  return random;
}

/* ----------------------------------------------------------------------
   integer RN between 1 and N inclusive
------------------------------------------------------------------------- */

int RandomFibLag::irandom(int n)
{
  int i = (int) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}
