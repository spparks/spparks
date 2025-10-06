/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS
RegionStyle(union,RegUnion)

#else

#ifndef SPK_REGION_UNION_H
#define SPK_REGION_UNION_H

#include "region.h"

namespace SPPARKS_NS {

class RegUnion : public Region {
 public:
  RegUnion(class SPPARKS *, int, char **);
  ~RegUnion();
  int match(double, double, double);

 private:
  int nregion;
  int *list;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Region union region ID does not exist

Self-explanatory.

*/
