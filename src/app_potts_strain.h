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
AppStyle(potts/strain,AppPottsStrain)

#else

#ifndef SPK_APP_POTTS_STRAIN_H
#define SPK_APP_POTTS_STRAIN_H

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsStrain : public AppPotts {
 public:
  AppPottsStrain(class SPPARKS *, int, char **);
  ~AppPottsStrain() {}
  void grow_app();
  void init_app();

  double site_propensity(int);
  void site_event(int, class RandomPark *);

  void *extract_app(char *);

 private:
  double *strain;
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

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
