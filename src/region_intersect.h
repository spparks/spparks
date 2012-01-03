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

#ifdef REGION_CLASS
RegionStyle(intersect,RegIntersect)

#else

#ifndef SPK_REGION_INTERSECT_H
#define SPK_REGION_INTERSECT_H

#include "region.h"

namespace SPPARKS_NS {

class RegIntersect : public Region {
 public:
  RegIntersect(class SPPARKS *, int, char **);
  ~RegIntersect();
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

E: Region intersect region ID does not exist

Self-explanatory.

*/
