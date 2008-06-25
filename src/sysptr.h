/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

// SysPtr class contains ptrs to master copy of
//   fundamental SPPARKS class ptrs stored in spparks.h
// top-level SPPARKS classes inherit from SysPtr to access spparks.h ptrs
// these variables are auto-initialized by SysPtr class constructor
// *& variables are really pointers to the pointers in spparks.h
// & enables them to be accessed directly in any class, e.g. error->all()

#ifndef SYSPTR_H
#define SYSPTR_H

#include "mpi.h"
#include "stdio.h"
#include "spparks.h"

namespace SPPARKS_NS {

class SysPtr {
public:
  explicit SysPtr(SPPARKS *ptr) :
    spk(ptr),
    universe(ptr->universe),
    input(ptr->input),
    memory(ptr->memory),
    error(ptr->error),
    app(ptr->app),
    solve(ptr->solve),
    sweep(ptr->sweep),
    timer(ptr->timer),
    output(ptr->output),
    world(ptr->world),
    infile(ptr->infile),
    screen(ptr->screen),
    logfile(ptr->logfile) {}
  virtual ~SysPtr() {}
  
 protected:
  SPPARKS *spk;

  Universe *&universe;      // universe of processors
  Input *&input;            // input script processing
  Memory *&memory;          // memory allocation functions
  Error *&error;            // error handling

  App *&app;                // application
  Solve *&solve;            // solver
  Sweep *&sweep;            // sweep
  Timer *&timer;            // timer
  Output *&output;          // output

  MPI_Comm &world;          // communicator for my world of procs
  FILE *&infile;            // infile for my world
  FILE *&screen;            // screen output for my world
  FILE *&logfile;           // logfile for my world
};

}

#endif


