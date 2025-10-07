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
AppStyle(potts/neigh,AppPottsNeigh)

#else

#ifndef SPK_APP_POTTS_NEIGH_H
#define SPK_APP_POTTS_NEIGH_H

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsNeigh : public AppPotts {
 public:
  AppPottsNeigh(class SPPARKS *, int, char **);
  ~AppPottsNeigh() {}
  void init_app();

  void site_event_rejection(int, class RandomPark *);
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
