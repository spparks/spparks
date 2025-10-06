/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(ising/single,AppIsingSingle)

#else

#ifndef SPK_APP_ISING_SINGLE_H
#define SPK_APP_ISING_SINGLE_H

#include "app_ising.h"

namespace SPPARKS_NS {

class AppIsingSingle : public AppIsing {
 public:
  AppIsingSingle(class SPPARKS *, int, char **);
  ~AppIsingSingle() {}

  void site_event_rejection(int, class RandomPark *);
};

}

#endif
#endif
