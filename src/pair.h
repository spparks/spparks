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

#ifndef SPK_PAIR_H
#define SPK_PAIR_H

#include "pointers.h"

namespace SPPARKS_NS {

class Pair : protected Pointers {
 public:
  int ntypes;
  double cutoff;

  Pair(class SPPARKS *);
  virtual ~Pair() {}
  void init();

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;
  virtual void init_style() {}
  virtual double init_one(int, int) {return 0.0;}
  virtual double energy(int, int, int *, double **, int *) = 0;

 protected:
  int allocated;                       // 0/1 = whether arrays are allocated
  int **setflag;
  double **cutsq;
  int mix_flag;

  double mix_energy(double, double, double, double);
  double mix_distance(double, double);
};

}

#endif

/* ERROR/WARNING messages:

E: All pair coeffs are not set

Self-explanatory.

*/
