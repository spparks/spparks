/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_SPPARKS_H
#define SPK_SPPARKS_H

#include "mpi.h"
#include "stdio.h"

namespace SPPARKS_NS {

class SPPARKS {
 public:
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling

  class App *app;                // application
  class Solve *solve;            // solver
  class Domain *domain;          // domain
  class Potential *potential;    // interatomic potentials
  class RanMars *ranmaster;      // master random number generator
  class Output *output;          // output
  class Timer *timer;            // timer

  MPI_Comm world;          // communicator for my world of procs
  FILE *infile;            // infile for my world
  FILE *screen;            // screen output for my world
  FILE *logfile;           // logfile for my world

  SPPARKS(int, char **, MPI_Comm);
  ~SPPARKS();
  void create();
  void destroy();

  void print_styles();
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch SPPARKS.

E: Processor partitions are inconsistent

The total number of processors in all partitions must match the number
of processors SPPARKS is running on.

E: Must use -in switch with multiple partitions

A multi-partition simulation cannot read the input script from stdin.
The -in command-line option must be used to specify a file.

E: Cannot open universe screen file

For a multi-partition run, the master screen file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open log.spparks

Self-explanatory.

E: Cannot open universe log file

For a multi-partition run, the master log file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open input script %s

Self-explanatory.

E: Cannot open screen file

The screen file specified as a command-line argument cannot be
opened.  Check that the directory you are running in allows for files
to be created.

E: Cannot open logfile

Self-explanatory.

E: Smallint setting in spktype.h is invalid

UNDOCUMENTED

E: Tagint setting in spktype.h is invalid

UNDOCUMENTED

E: Bigint setting in spktype.h is invalid

UNDOCUMENTED

E: MPI_SPK_TAGINT and tagint in spktype.h are not compatible

UNDOCUMENTED

E: MPI_SPK_BIGINT and bigint in spktype.h are not compatible

UNDOCUMENTED

*/
