/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "pointers.h"

namespace SPPARKS_NS {

class Output : protected Pointers {
 public:
  Output(class SPPARKS *);
  ~Output();
  void init(double);
  double setup(double);
  double compute(double, int);
  void set_stats(int, char **);
  void add_dump(int, char **);
  void dump_one(int, char **);
  void dump_modify(int, char **);
  void undump(int, char **);
  void add_diag(class Diag *);

 private:
  int me,nprocs;

  double stats_time,stats_delta;     // stats info
  double stats_scale,stats_delay;
  int stats_logfreq,stats_nrepeat;

  int ndump;                         // list of dumps
  int max_dump;
  class Dump **dumplist;

  int ndiag;                         // list of diagnostics
  class Diag **diaglist;

  void stats(int);
  void stats_header();
  double next_time(double, int, double, int, double, double);
};

}

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Reuse of dump ID

UNDOCUMENTED

E: Invalid dump style

UNDOCUMENTED

E: Could not find dump ID in dump_one command

Self-explanatory.

E: Cannot use dump_one for first snapshot in dump file

Self-explanatory.

E: Could not find dump ID in dump_modify command

Self-explanatory.

E: Could not find dump ID in undump command

Self-explanatory.

*/
