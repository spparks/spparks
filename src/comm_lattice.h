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
  void init(class SweepLattice *, const int, const int, int *);
  void all();
  void sector(int);
  void reverse_sector(int) {}

 private:
  int me,nprocs;
  int delghost,dellocal; // thickness of ghost and local communication layers 

  struct Swap {
    int nsend,nrecv;               // number of messages to send/recv
    int *sproc;                    // proc for each send message
    int *scount;                   // size of each send message in sites
    int *smax;                     // max size of each send message in sites
    int **sindex;                  // list of my lattice indices for each send
    int *sibuf;                    // biggest int send message
    double *sdbuf;                 // biggest double send message
    int *rproc;                    // proc for each recv message
    int *rcount;                   // size of each recv message in sites
    int *rmax;                     // max size of each recv message in sites
    int **rindex;                  // list of my lattice indices for each recv
    int **ribuf;                   // each int recv message
    double **rdbuf;                // each double recv message
    MPI_Request *request;          // MPI datums for each recv message
    MPI_Status *status;
  };

  Swap *allswap;
  Swap **sectorswap;
  int nsector;

  struct Ghost {
    int id_global,index_local,proc;
  };

  int sitecustom;
  int ninteger,ndouble;
  int *lattice;
  int **iarray;
  double **darray;

  Swap *create_swap(int, int *, int, int *, int *, int **);
  void free_swap(Swap *);
  void perform_swap_lattice(Swap *);
  void perform_swap_int(Swap *);
  void perform_swap_double(Swap *);
  void perform_swap_general(Swap *);
};

}

#endif
