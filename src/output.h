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

namespace SPPARKS_NS {

class Output : protected Pointers {
 public:
  Output(class SPPARKS *);
  ~Output();
  void init(double);
  double setup(double);
  void set_stats(int, char **);
  void set_dump(int, char **);
  void add_diag(class Diag *);
  double compute(double, int);
  void dump_one(double);

 private:
  int me,nprocs;

  double stats_time,stats_delta,stats_scale,stats_delay;
  int stats_logfreq,stats_nrepeat;

  double dump_time,dump_delta,dump_scale,dump_delay;
  int dump_logfreq,dump_nrepeat;
  int idump;

  int ndiags;
  class Diag **diaglist;

  FILE *fp;                  // dump file pointer
  int size_one;

  class AppLattice *applattice;
  class AppLattice2d *applattice2d;
  class AppLattice3d *applattice3d;

  int nglobal,nlocal,nx_local,ny_local,nz_local;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;
  int mask_flag;

  int *vtype;                // type of each vector (INT, DOUBLE)
  int *vindex;               // index into int,double packs
  char **vformat;            // format string for each vector element

  double *buf;
  int *mask;
  int maxbuf;

  double next_time(double, int, double, int, double, double);
  void stats(int);
  void stats_header();
  void dump_header();
  void dump(double);
  void write_data(int, double *);

  typedef void (Output::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_lattice(int);
  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_energy(int);
  void pack_propensity(int);
  void pack_integer(int);
  void pack_double(int);

  void maskzeroenergy();
};

}

#endif
