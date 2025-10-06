/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(propensity,DiagPropensity)

#else

#ifndef SPK_DIAG_PROPENSITY_H
#define SPK_DIAG_PROPENSITY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagPropensity : public Diag {

 public:
  DiagPropensity(class SPPARKS *, int, char **);
  ~DiagPropensity() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppLattice *applattice;
  int nlocal;
  double propensity;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

E: Diag propensity requires KMC solve be performed

Only KMC solvers compute a propensity for sites and the system.

*/
