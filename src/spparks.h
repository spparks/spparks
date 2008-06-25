/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SPPARKS_H
#define SPPARKS_H

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
  class Sweep *sweep;            // sweep
  class Timer *timer;            // timer
  class Output *output;          // output

  MPI_Comm world;          // communicator for my world of procs
  FILE *infile;            // infile for my world
  FILE *screen;            // screen output for my world
  FILE *logfile;           // logfile for my world

  SPPARKS(int, char **, MPI_Comm);
  ~SPPARKS();
  void create();
  void init();
  void destroy();
};

}

#endif
