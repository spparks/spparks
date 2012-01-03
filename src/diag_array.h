/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
 
   This class added by Eric Homer (copied after diag_energy.h).

   See the README file in the top-level SPPARKS directory.
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
