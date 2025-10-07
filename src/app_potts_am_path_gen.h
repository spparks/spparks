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
AppStyle(potts/am/path/gen,AppPottsAmPathGen)

#else

#ifndef SPK_APP_POTTS_AM_PATH_GEN
#define SPK_APP_POTTS_AM_PATH_GEN

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

class AppPottsAmPathGen : public PottsAmPathParser {

   public:
      AppPottsAmPathGen(class SPPARKS *spk, int narg, char **arg);
      virtual ~AppPottsAmPathGen();
      virtual void grow_app() {}
      virtual void init_app();
      void input_app(char *, int , char **);
      void app_update(double);

   private:
      std::string path_filename="";
      int num_layers=-1, zstart=0;
      int melt_depth=0, width_haz=0, depth_haz=0;
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
