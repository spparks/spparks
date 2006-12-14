/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef UNIVERSE_H
#define UNIVERSE_H

#include "mpi.h"
#include "stdio.h"
#include "sysptr.h"

namespace SPPARKS {

class Universe : protected SysPtr {
 public:
  char *version;          // SPPARKS version string = date

  MPI_Comm uworld;        // communicator for entire universe
  int me,nprocs;          // my place in universe

  FILE *uscreen;          // universe screen output
  FILE *ulogfile;         // universe logfile

  int nworlds;            // # of worlds in universe
  int iworld;             // which world I am in
  int *procs_per_world;   // # of procs in each world
  int *root_proc;         // root proc in each world

  explicit Universe(class SPK *, MPI_Comm);
  ~Universe();
  void add_world(char *);
  int consistent();

 private:
  Universe(); // Not a sane operation.
  Universe(const Universe&); // Not a sane operation.
  Universe& operator=(const Universe&); // Not a sane operation.
};

}

#endif
