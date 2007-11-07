/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_LATTICE_H
#define COMM_LATTICE_H

#include "mpi.h"
#include "sysptr.h"

namespace SPPARKS {

class CommLattice : protected SysPtr {
 public:
  CommLattice(class SPK *);
  ~CommLattice();
  void init(class SweepLattice *, const int, const int);
  void all(int *);
  void sector(int *, const int);
  void reverse_sector(int *, const int) {}

 private:
  int me,nprocs;
  int dimension;
  int nlocal;
  int procwest,proceast,procsouth,procnorth,procdown,procup;
  int delghost,dellocal; // thickness of ghost and local communication layers 

  struct Swap {
    int nsend,nrecv;               // number of messages to send/recv
    int *sproc;                    // proc for each send message
    int *scount;                   // size of each send message
    int **sindex;                  // list of my lattice indices for each send
    int *sbuf;                     // length of biggest send message
    int *rproc;                    // proc for each recv message
    int *rcount;                   // size of each receive message
    int **rindex;                  // list of my lattice indices for each recv
    int **rbuf;                    // incoming buf for each recv message
    MPI_Request *request;          // MPI datums for each recv message
    MPI_Status *status;
  };

  Swap *allswap;
  Swap **sectorswap;
  int nsector;

  struct Ghost {
    int id,proc,index;
  };

  Swap *create_swap(int, int *, int, int *, int *, int **);
  void free_swap(Swap *);
  void perform_swap(Swap *, int *);
};

}

#endif

