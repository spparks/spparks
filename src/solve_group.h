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

#ifdef SOLVE_CLASS
SolveStyle(group,SolveGroup)

#else

#ifndef SPK_SOLVE_GROUP_H
#define SPK_SOLVE_GROUP_H

#include "solve.h"

namespace SPPARKS_NS {

class SolveGroup : public Solve {
 public:
  SolveGroup(class SPPARKS *, int, char **);
  ~SolveGroup();
  SolveGroup *clone();

  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  class Groups *groups;
  int nevents;

  double *p;                     // local copy of propensities
  double lo,hi;
  int ngroups;

  int nroundlo,nroundhi;         // info on propensities reset to lo/hi
  double lomax,himax;

  void round_check();
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

E: Invalid probability bounds for solve_style group

Self-explanatory.

W: %d propensities were reset to lo value, max lo = %g

UNDOCUMENTED

W: %d propensities were reset to hi value, max hi = %g

UNDOCUMENTED

*/
