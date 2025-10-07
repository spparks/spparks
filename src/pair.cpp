/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};
enum{R,RSQ,BMP};

/* ---------------------------------------------------------------------- */

Pair::Pair(SPPARKS *spk) : Pointers(spk)
{
  mix_flag = GEOMETRIC;
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
  int i,j;

  if (!allocated) error->all(FLERR,"All pair coeffs are not set");

  // I,I coeffs must be set
  // init_one() will check if I,J is set explicitly or inferred by mixing

  for (i = 1; i <= ntypes; i++)
    if (setflag[i][i] == 0) error->all(FLERR,"All pair coeffs are not set");

  // style-specific initialization

  init_style();

  // call init_one() for each I,J
  // set cutsq for each I,J, used to neighbor
  // cutforce = max of all I,J cutoffs

  double cut;
  cutoff = 0.0;
  for (i = 1; i <= ntypes; i++)
    for (j = i; j <= ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
      cutoff = MAX(cutoff,cut);
    }
}

/* ----------------------------------------------------------------------
   mixing of pair potential prefactors (epsilon)
------------------------------------------------------------------------- */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == ARITHMETIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == SIXTHPOWER)
    value = 2.0 * sqrt(eps1*eps2) *
      pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) * pow(sig2,6.0));
  return value;
}

/* ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
------------------------------------------------------------------------- */

double Pair::mix_distance(double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(sig1*sig2);
  else if (mix_flag == ARITHMETIC)
    value = 0.5 * (sig1+sig2);
  else if (mix_flag == SIXTHPOWER)
    value = pow((0.5 * (pow(sig1,6.0) + pow(sig2,6.0))),1.0/6.0);
  return value;
}
