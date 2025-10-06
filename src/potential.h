/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_POTENTIAL_H
#define SPK_POTENTIAL_H

#include "pointers.h"

namespace SPPARKS_NS {

class Potential : protected Pointers {
 public:
  char *pair_style;
  class Pair *pair;

  Potential(class SPPARKS *);
  ~Potential();
  void init();
  void create_pair(const char *);
  void bounds(char *, int, int &, int &);

 private:
  class Pair *new_pair(const char *);
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid pair style

Self-explanatory.

*/
