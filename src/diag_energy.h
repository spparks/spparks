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

#ifdef DIAG_CLASS
DiagStyle(energy,DiagEnergy)

#else

#ifndef SPK_DIAG_ENERGY_H
#define SPK_DIAG_ENERGY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagEnergy : public Diag {
 public:
  DiagEnergy(class SPPARKS *, int, char **);
  ~DiagEnergy() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appofflattice;
  double energy;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

*/
