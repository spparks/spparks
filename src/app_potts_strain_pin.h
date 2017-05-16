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

#ifdef APP_CLASS
AppStyle(potts/strain/pin,AppPottsStrainPin)

#else

#ifndef SPK_APP_POTTS_STRAIN_PIN_H
#define SPK_APP_POTTS_STRAIN_PIN_H

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsStrainPin : public AppPotts {
 public:
  AppPottsStrainPin(class SPPARKS *, int, char **);
  ~AppPottsStrainPin() {}
  void grow_app();
  void init_app();
  void input_app(char *, int, char **);

  double site_propensity(int);
  void site_event(int, class RandomPark *);
  void push_new_site(int, int*, int, std::stack<int>*);

  void *extract_app(char *);

 private:
  double *strain;
  double pfraction;
  int multi,nthresh;

  void pin_create();
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
