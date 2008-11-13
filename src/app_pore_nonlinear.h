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

#ifndef APP_PORE_NONLINEAR_H
#define APP_PORE_NONLINEAR_H

#include "app_diffusion_nonlinear.h"

namespace SPPARKS_NS {

class AppPoreNonLinear : public AppDiffusionNonLinear {
 public:
  AppPoreNonLinear(class SPPARKS *, int, char **);
};

}

#endif
