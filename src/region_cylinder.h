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
