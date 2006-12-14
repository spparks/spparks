/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SPK_H
#define SPK_H

#include "mpi.h"
#include "stdio.h"

namespace SPPARKS {

class SPK {
 public:
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling

  class App *app;                // application
  class Solve *solve;            // solver
  class Timer *timer;            // timer

  MPI_Comm world;          // communicator for my world of procs
  FILE *infile;            // infile for my world
  FILE *screen;            // screen output for my world
  FILE *logfile;           // logfile for my world

  SPK(int, char **, MPI_Comm);
  ~SPK();
  void create();
  void init();
  void destroy();
};

}

#endif
