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
DiagStyle(erbium,DiagErbium)

#else

#ifndef SPK_DIAG_ERBIUM_H
#define SPK_DIAG_ERBIUM_H

#include "diag.h"

namespace SPPARKS_NS {

class DiagErbium : public Diag {
 public:
  DiagErbium(class SPPARKS *, int, char **);
  ~DiagErbium();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppErbium *apperbium;
  int nlist;
  char **list;
  int *which,*index,*ivector;
  int siteflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag_style erbium requires app_style erbium

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Invalid value setting in diag_style erbium

UNDOCUMENTED

*/
