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
AppStyle(potts/grad,AppPottsGrad)

#else

#ifndef SPK_APP_POTTSGRAD_H
#define SPK_APP_POTTSGRAD_H


#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsGrad : public AppPotts {
 public:
   AppPottsGrad(class SPPARKS *, int, char **);
   virtual ~AppPottsGrad();
   virtual void init_app();
   virtual void grow_app();
   double compute_temperature(double x, double y, double z) const;
   double compute_site_temperature(int site_i) const;
   double compute_mobility(double temperature_site_i, double max_T) const;
   double mobility_grad(int i);

   virtual void site_event_rejection(int, class RandomPark *);

 protected:

   // temperature @ centroid
   double T0;
   // temperature gradient @ centroid
   double grad_x, grad_y, grad_z;
   // maximum temperature on domain
   double max_T;
   // maximum mobility
   double max_M;
   // Site temperatures; T0 + grad_tx * dx + grad_ty * dy + grad_tz * dz
   double *T;
   // mobility m=m0*exp(-q/kT)
   double m0;
   // mobility values for each site
   double *M;
   // activation energy q
   double activation_energy;
   // length scale conversion factor
   double convert;
   // determines if temp gradient or mobility gradient is used
   char *gradient_choice;

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
