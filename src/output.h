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

#ifndef SPK_OUTPUT_H
#define SPK_OUTPUT_H

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
  void dump_one(int, char **, double);
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

#endif
