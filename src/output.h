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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "pointers.h"
#include "diag.h"

namespace SPPARKS_NS {

class Output : protected Pointers {
 public:
  Output(class SPPARKS *);
  ~Output();

  void init(double);
  void set_stats(int, char **);
  void set_dump(int, char **);
  void stats();
  void stats_header();
  void dump_header();
  void compute(double, int);
  void add_diag(Diag *);
  int ndiags;
  Diag** diaglist;

 private:
  double stats_time,stats_delta,stats_scale,stats_t0;
  double dump_time,dump_delta;
  int me,nprocs;
  int stats_nrepeat,stats_irepeat,stats_ilogfreq;
};

}

#endif
