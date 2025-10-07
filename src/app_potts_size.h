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
// clang-format off
AppStyle(potts/size,AppPottsSize)
// clang-format on
#else

#ifndef SPK_APP_POTTS_SIZE_H
#define SPK_APP_POTTS_SIZE_H

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsSize : public AppPotts {
public:
  AppPottsSize(class SPPARKS *, int, char **);
  ~AppPottsSize() {}
  void grow_app();
  void init_app();

  double site_propensity(int);
  void site_event(int, class RandomPark *);

  void *extract_app(char *);

private:
  double *size;
};

} // namespace SPPARKS_NS

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
