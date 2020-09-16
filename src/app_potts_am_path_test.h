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
AppStyle(potts/am/path/test,AppPottsAmPathTest)

#else

#ifndef SPK_APP_POTTS_AM_PATH_TEST
#define SPK_APP_POTTS_AM_PATH_TEST

#include <vector>
#include <stdlib.h>
#include <tuple>
#include <map>
#include "random_park.h"
#include "app_potts.h"
#include "am_raster.h"
#include "potts_am_path_parser.h"

using std::vector;
using std::tuple;

namespace SPPARKS_NS {

class AppPottsAmPathTest : public PottsAmPathParser {

   typedef tuple<double,double,double,double> ComputationalVolume;

   public:
      AppPottsAmPathTest(class SPPARKS *, int, char **);
      virtual ~AppPottsAmPathTest();
      virtual void grow_app() {}
      virtual void init_app();
      void input_app(char *, int , char **);
      void app_update(double);
      void print_path();
      tuple<vector<double>,vector<ComputationalVolume>>
      get_layer_computational_volumes(const Pass& p, START s, int offset_x, int offset_y, int width_haz) const;

   private:
      std::string path_filename="";
      int num_layers=-1;
      int melt_depth=0, depth_haz=0;
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
