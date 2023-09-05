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
AppStyle(potts/am/bezier,AppPottsAmBezier)

#else

#ifndef SPK_APP_POTTS_AM_BEZIER
#define SPK_APP_POTTS_AM_BEZIER

#include <stdlib.h>
#include <array>
#include "app_potts.h"
#include "am_raster.h"
#include "potts_am_path_parser.h"
#include "am_bezier.h"
#include "fourth_order_bezier.h"

namespace SPPARKS_NS {

using std::array;

class AppPottsAmBezier : public PottsAmPathParser {

public:
   AppPottsAmBezier(class SPPARKS *, int, char **);
   virtual void grow_app();
   virtual void init_app();
   virtual void site_event_rejection(int, RandomPark *);
   void input_app(char *, int , char **);
   double compute_mobility(int) const;
   void app_update(double);


  // private functions
private:
   void parse_am_bezier(int narg, char** arg);
   array<size_t,3> get_ijk(size_t index, size_t nx, size_t ny) const;
   std::vector<double> compute_dLookup_table
      (int num_points, size_t i0, size_t j0, size_t k0, 
       const std::vector<double>& dx, const std::vector<double>& dz) const;
   // private variables
private:
   double pool_width, pool_depth, haz; 
   array<double,5> xt, yt, zs;
   double beta_y, beta_z;
   double *distance;
   AM_Bezier bezier;
   RASTER::rectangular_range HAZ;
   std::vector<double> d_table;
   std::vector<double> dim_d_table={0,0,0};
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
