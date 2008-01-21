/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef OUTPUT_H
#define OUTPUT_H

#include "sysptr.h"
#include "diag.h"
namespace SPPARKS {

class Output : protected SysPtr {
 public:
  explicit Output(class SPK *);
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
  double stats_time,stats_delta,stats_scale;
  double dump_time,dump_delta;
  int stats_nrepeat, stats_irepeat, stats_ilogfreq,stats_t0;
};

}

#endif
