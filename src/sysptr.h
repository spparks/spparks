/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   In the following description, classname System refers to the
   top-level class, which in this case is called SPK.
   Parent class for all classes in the code, except System.
   Contains references to member data of an System object.
   System contains a single instance of each of the major classes.
   This construction prevents multiple System instances from 
      stomping on each other while avoiding the need to have
      each object in the code explicitly identify which System 
      object they belong to.
------------------------------------------------------------------------- */

#ifndef SYSPTR_H
#define SYSPTR_H

#include "mpi.h"
#include "spk.h"

namespace SPPARKS {

class SysPtr {
public:
  explicit SysPtr(SPK *sys) :
    spk(sys),
    universe(sys->universe),
    input(sys->input),
    memory(sys->memory),
    error(sys->error),
    app(sys->app),
    solve(sys->solve),
    sweep(sys->sweep),
    timer(sys->timer),
    output(sys->output),
    world(sys->world),
    infile(sys->infile),
    screen(sys->screen),
    logfile(sys->logfile) {}
  virtual ~SysPtr() {}

protected:
  SPK *spk;

  class Universe *&universe;      // universe of processors
  class Input *&input;            // input script processing
  class Memory *&memory;          // memory allocation functions
  class Error *&error;            // error handling

  class App *&app;                // application
  class Solve *&solve;            // solver
  class Sweep *&sweep;            // sweep
  class Timer *&timer;            // timer
  class Output *&output;          // output

  MPI_Comm &world;          // communicator for my world of procs
  FILE *&infile;            // infile for my world
  FILE *&screen;            // screen output for my world
  FILE *&logfile;           // logfile for my world
};

}

#endif


