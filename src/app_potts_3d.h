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

#ifndef APP_POTTS_3D_H
#define APP_POTTS_3D_H

#include "app_lattice3d.h"

namespace SPPARKS_NS {

class AppPotts3d : public AppLattice3d {
 public:
  AppPotts3d(class SPPARKS *, int&, char **);
  ~AppPotts3d();

 protected:
  int nspins;
  char * spinfile;
};

}

#endif
