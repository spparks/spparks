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
