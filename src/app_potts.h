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
AppStyle(potts,AppPotts)

#else

#ifndef SPK_APP_POTTS_H
#define SPK_APP_POTTS_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppPotts : public AppLattice {
 public:
  AppPotts(class SPPARKS *, int, char **);
  virtual ~AppPotts();
  virtual void grow_app();
  virtual void init_app();

  virtual double site_energy(int);
  virtual void site_event_rejection(int, class RandomPark *);
  virtual double site_propensity(int);
  virtual void site_event(int, class RandomPark *);

 protected:
  int nspins;
  int *spin;
  int *sites,*unique;
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
