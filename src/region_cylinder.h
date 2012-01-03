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
RegionStyle(cylinder,RegCylinder)

#else

#ifndef SPK_REGION_CYLINDER_H
#define SPK_REGION_CYLINDER_H

#include "region.h"

namespace SPPARKS_NS {

class RegCylinder : public Region {
 public:
  RegCylinder(class SPPARKS *, int, char **);
  int match(double, double, double);

 private:
  char axis;
  double c1,c2;
  double radius;
  double lo,hi;
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

E: Cannot use region INF or EDGE when box does not exist

Can only define a region with these parameters after a simulation
box has been defined.

*/
