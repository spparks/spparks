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

#ifndef SPK_DUMP_H
#define SPK_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace SPPARKS_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  int idump;                 // counter for output snapshots
  double next_time,delta;    // params governing output times
  double scale,delay;
  int logfreq,nrepeat;

  Dump(class SPPARKS *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write(double);
  void modify_params(int, char **);

 protected:
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int padflag;               // timestep padding in filename

  double *buf;               // memory for site quantities
  int maxbuf;                // size of buf
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one site

  int latticeflag;           // 1 for on-lattice, 0 for off-lattice app
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;

  virtual void init_style() = 0;
  void openfile();
  virtual int modify_param(int, char **) {return 0;}
  virtual void write_header(int, double) = 0;
  virtual int count() = 0;
  virtual void pack() = 0;
  virtual void write_data(int, double *) = 0;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Dump command can only be used for spatial applications

Self-explanatory.

E: Cannot open gzipped file

Self-explantory.

E: Cannot open dump file

Self-explanatory.

*/
