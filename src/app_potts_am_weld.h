/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts/am/weld,AppPottsAmWeld)

#else

#ifndef SPK_APP_POTTS_AM_WELD
#define SPK_APP_POTTS_AM_WELD

#include <vector>
#include "random_park.h"
#include <stdlib.h>
#include <map>
#include "am_raster.h"
#include "potts_am_path_parser.h"

using RASTER::pool_shape::ShapeType;

namespace SPPARKS_NS {

class AppPottsAmWeld : public PottsAmPathParser {

 public:
  AppPottsAmWeld(class SPPARKS *, int, char **);
  virtual ~AppPottsAmWeld() { }
  virtual void grow_app();
  virtual void init_app();
  virtual void site_event_rejection(int, RandomPark *);
  void input_app(char *, int , char **);
  double compute_mobility(int);
  void app_update(double);

 private:
   double alpha, beta;
   double haz;
   double *distance;
   SPPARKS_NS::RandomPark random_park;
   double simulation_time;

   // Pool shape parameters
   ShapeType shape_type;
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
