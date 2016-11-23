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
  int firstflag;             // 1 if no dump yet written, 0 if not first
  double next_time,delta;    // params governing output times
  double scale,delay;
  int logfreq,nrepeat;

  // static variable across all Dump objects

  static Dump *dumpptr;         // holds a ptr to Dump currently being used

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
                             // else # of procs writing files
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  char *multiname;           // filename with % converted to cluster ID
  MPI_Comm clustercomm;      // MPI communicator within my cluster of procs

  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int padflag;               // timestep padding in filename
  int singlefile_opened;     // 1 = one big file, already opened, else 0

  int sort_flag;             // 1 if sorted output
  int sortcol;               // 0 to sort on ID, 1-N on columns
  int sortcolm1;             // sortcol - 1
  int sortorder;             // ASCEND or DESCEND

  int maxbuf;                // size of buf
  double *buf;               // memory for site quantities

  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one site
  int nme;                   // # of sites in this dump from me
  bigint ntotal;             // total # of per-site lines in snapshot

  int maxids;                // size of ids
  int maxsort;               // size of bufsort, idsort, index
  int maxproc;               // size of proclist
  int reorderflag;           // 1 if OK to reorder instead of sort
  int ntotal_reorder;        // # of sites that must be in snapshot
  int nme_reorder;           // # of sites I must own in snapshot
  tagint idlo;               // lowest ID I own when reordering

  tagint *ids;               // list of site IDs, if sorting on IDs
  double *bufsort;
  tagint *idsort;
  int *index,*proclist;

  class Irregular *irregular;

  int latticeflag;           // 1 for on-lattice, 0 for off-lattice app
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;

  virtual void init_style() = 0;
  void openfile();
  virtual int modify_param(int, char **) {return 0;}
  virtual void write_header(bigint, double) = 0;
  virtual void write_footer() {}
  virtual int count() = 0;
  virtual void pack(tagint *) = 0;
  virtual void write_data(int, double *) = 0;

  void sort();
  static int idcompare(const void *, const void *);
  static int bufcompare(const void *, const void *);
  static int bufcompare_reverse(const void *, const void *);
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
