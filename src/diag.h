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

#ifndef SPK_DIAG_H
#define SPK_DIAG_H

#include "pointers.h"

namespace SPPARKS_NS {

class Diag : protected Pointers {
 public:
  char *style;
  int stats_flag;                   // 1 if stats drives output, 0 if not
  double next_time,delta;           // output params for stats_flag = 0
  double scale,delay;
  int logfreq,nrepeat;

  Diag(class SPPARKS *, int, char **);
  virtual ~Diag();

  // pure virtual functions, must be defined in child class
  
  virtual void init() = 0;
  virtual void compute() = 0;

  // virtual functions, may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};

 protected:
  int me,nprocs;
  int iarg_child;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

*/
