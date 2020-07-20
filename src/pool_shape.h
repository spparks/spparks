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

#ifndef SPK_APP_POTTSWELD_POOL_SHAPE_H
#define SPK_APP_POTTSWELD_POOL_SHAPE_H

#include<vector>

using std::vector;

namespace weld {

namespace pool_shape {
   enum ShapeType {undefined=-1, ellipse=0, teardrop=1};

   class PoolShape {
      public:
         virtual double distance(const double *XYZ) const = 0;
         PoolShape() {}
         virtual ~PoolShape() {}
   };

}

}

#endif
