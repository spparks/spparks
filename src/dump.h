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

#include "stdio.h"
#include "pointers.h"

namespace SPPARKS_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  int idump;                 // counter for output snapshots
  double next_time,delta;    // params governing output times
  double scale,delay;
  int logfreq,nrepeat;

  Dump(class SPPARKS *, int, char **);
  ~Dump();
  void init();
  void modify_params(int, char **);
  void write(double);

 private:
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc
  int flush_flag;            // 0 if no flush, 1 if flush every dump

  double *buf;               // memory for atom quantities
  int maxbuf;                // size of buf
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one atom

  int *vtype;                // type of each vector (INT, DOUBLE)
  int *vindex;               // index into int,double packs
  char **vformat;            // format string for each vector element

  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshold on for each nthresh
  int *thresh_op;            // threshold operation for each nthresh
  double *thresh_value;      // threshold value for each nthresh
  int *thresh_index;         // N index for iN and dN thresholds

  char *columns;             // text describing columns of dump output

  int *choose;               // lists of sites chosen for output
  double *dchoose;
  int maxlocal;              // size of choose arrays

  int latticeflag;           // 1 for on-lattice, 0 for off-lattice app
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  int nglobal,nlocal,nx_local,ny_local,nz_local;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;
  int mask_flag;

  int count();
  void pack();
  void write_header(int, double);
  void write_data(int, double *);
  void openfile();

  typedef void (Dump::*FnPtrHeader)(int, double);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int, double);
  void header_text(int, double);

  typedef void (Dump::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_text(int, double *);

  typedef void (Dump::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_site(int);
  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_energy(int);
  void pack_propensity(int);
  void pack_iarray(int);
  void pack_darray(int);
};

}
