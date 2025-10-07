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
DiagStyle(array,DiagArray)

#else

#ifndef SPK_DIAG_ARRAY_H
#define SPK_DIAG_ARRAY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagArray : public Diag {
 public:
  DiagArray(class SPPARKS *, int, char **);
  ~DiagArray();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appofflattice;
  double *vals;
  int *index;
  bool *index_double;       // false for int, true for double
  int *diag_method;
  int nvals;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

E: Invalid diag_style command

UNDOCUMENTED

*/
