/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_EPROF3D_H
#define DIAG_EPROF3D_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS {

class DiagEprof3d : public Diag {
  friend class SweepLattice3d;

 public:
  DiagEprof3d(class SPK *, int, char **);
  virtual ~DiagEprof3d();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  // Functions and Data for Eprof Analysis
  void write_header();
  void write_prof(double);

  double*** lat3d;
  int ndata;
  double *prof;
  double *count;
  int *ixtable, *iytable, *iztable;

  FILE *fp;
  class AppLattice3d *applattice3d;

  int nx_global,ny_global,nz_global,nx_procs,ny_procs,nz_procs,delghost;
  int nx_local,ny_local,nz_local,nx_offset,ny_offset,nz_offset;
  int nxlo,nylo,nzlo,nxhi,nyhi,nzhi;

  int prof_index;

  enum ProfStyles {STANDARD,BOUNDARY};
  int prof_style;
  int iboundary;
  int nbound;
  double eav,eb1,eb2;

};

}

#endif
