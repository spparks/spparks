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

#ifndef DIAG_EPROF3D_H
#define DIAG_EPROF3D_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagEprof3d : public Diag {
  friend class SweepLattice3d;

 public:
  DiagEprof3d(class SPPARKS *, int, char **);
  virtual ~DiagEprof3d();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  // Functions and Data for Eprof Analysis
  void write_header();
  void write_prof(double);

  int ndata;
  double *prof;
  double *count;
  int *ixtable, *iytable, *iztable;

  FILE *fp;
  class AppLattice3d *applattice3d;

  int nx_global,ny_global,nz_global,nx_procs,ny_procs,nz_procs;
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
