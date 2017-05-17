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
AppStyle(potts/weld,AppPottsWeld)

#else

#ifndef SPK_APP_POTTSWELD_H
#define SPK_APP_POTTSWELD_H

#include <vector>
#include "random_park.h"
#include "app_potts.h"
#include "pool_shape.h"

using std::vector;

namespace SPPARKS_NS {

class AppPottsWeld : public AppPotts {
 public:
   AppPottsWeld(class SPPARKS *, int, char **);
   virtual ~AppPottsWeld() { }
   virtual void init_app();
   virtual void grow_app();
   virtual void site_event_rejection(int, class RandomPark *);
   virtual void app_update(double dt);
   virtual void input_app(char *, int, char **);
   double compute_mobility(int site) const;
   void print_potts_weld_params() const;

 private:
   double yp;
   double alpha, beta;
   double velocity;
   double haz;
   double *distance;
   SPPARKS_NS::RandomPark random_park;
   double simulation_time;
   double pulse_amplitude, pulse_step_frequency;

   // Pool shape parameters
   weld::pool_shape::ShapeType shape_type;
   double width, length;
   vector<vector<double> > teardrop_control_points;

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
