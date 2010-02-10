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

#ifndef SPK_COMM_OFF_LATTICE_H
#define SPK_COMM_OFF_LATTICE_H

#include "mpi.h"
#include "pointers.h"

namespace SPPARKS_NS {

class CommOffLattice : protected Pointers {
 public:
  CommOffLattice(class SPPARKS *);
  ~CommOffLattice();
  void init(int);
  void all();
  void all_reverse();
  void sector(int);
  void reverse_sector(int);

 private:
  int me,nprocs;

  struct Swap {
    int nsend,nrecv;               // number of messages to send/recv
    int *sproc;                    // proc for each send message
    int *scount;                   // # of bins to send each proc
    int **sindex;                  // list of my bin indices for each send
    int **ssite;                   // list of site counts, one per send bin
    int *stotal;                   // total # of sites to send to each proc
    int *rproc;                    // proc for each recv message
    int *rcount;                   // # of bins to recv from each proc
    int **rindex;                  // list of my bin indices for each recv
    int **rsite;                   // list of site counts, one per recv bin
    int *rtotal;                   // total # of sites to recv from each proc
    int ncopy;                     // 1 if copy bins to self, 0 if not
    int ccount;                    // # of bins to copy
    int *cbinsrc;                  // list of my bin indices to copy from
    int *cbindest;                 // list of my bin indices to copy to
    MPI_Request *request;          // MPI datums for each recv message
    MPI_Status *status;
  };

  Swap *allswap;
  Swap **sectorswap;
  Swap **sectorreverseswap;
  int nsector;

  class AppOffLattice *appoff;

  int size_one;
  int smax,rmax;
  double *sbuf,*rbuf;

  int ninteger,ndouble;
  int site_only;                    // only 1 int, no doubles
  double xprd,yprd,zprd;

  int nchunk;
  int chunklo[8][5][3],chunkhi[8][5][3];

  Swap *create_swap_all();
  Swap *create_swap_sector(int);
  Swap *create_swap_sector_reverse(int);
  void free_swap(Swap *);

  void create_send_from_list(int, int **, Swap *);
  void create_recv_from_send(Swap *);

  void perform_swap(Swap *);
  void perform_swap_reverse(Swap *);

  int bin_sector_ghost(int, int, int, int, int);
  void setup_sector_chunks(int, int, int);
  void fill_chunk(int, int, int, int, int, int, int, int);
};

}

#endif
