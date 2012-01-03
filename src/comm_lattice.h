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

#ifndef SPK_COMM_LATTICE_H
#define SPK_COMM_LATTICE_H

#include "mpi.h"
#include "pointers.h"

namespace SPPARKS_NS {

class CommLattice : protected Pointers {
 public:
  CommLattice(class SPPARKS *); ~CommLattice();
  void init(int, int, int, int *);
  void all();
  void all_reverse();
  void sector(int);
  void reverse_sector(int);

 private:
  int me,nprocs;
  int delghost,delreverse;

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

  struct Site {
    tagint id_global;
    int index_local;
    int proc;
  };

  Swap *allswap;
  Swap *reverseswap;
  Swap **sectorswap;
  Swap **sectorreverseswap;
  int nsector;

  int ninteger,ndouble;
  int site_only;                     // only 1 int, no doubles
  int **iarray;
  double **darray;
  int *site;                         // simply points to iarray[0]

  Swap *create_swap_all();
  Swap *create_swap_all_reverse();
  Swap *create_swap_sector(int, int *);
  Swap *create_swap_sector_reverse(int, int *);
  void free_swap(Swap *);

  void create_send_from_list(int, Site *, Swap *);
  void create_send_from_recv(int, int, Site *, Swap *);
  void create_recv_from_send(int, int, Site *, Swap *);
  void create_recv_from_list(int, Site *, Swap *);

  void perform_swap_site(Swap *);
  void perform_swap_int(Swap *);
  void perform_swap_double(Swap *);
  void perform_swap_general(Swap *);
};

}

#endif

/* ERROR/WARNING messages:

E: Site-site interaction was not found

Internal SPPARKS error.

*/
