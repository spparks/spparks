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
    int *scount;                   // size of each send message in bins
    int **sindex;                  // list of my bin indices for each send
    int *rproc;                    // proc for each recv message
    int *rcount;                   // size of each recv message in bins
    int **rindex;                  // list of my bin indices for each recv
    int ncopy;                     // 1 if copy bins to self, 0 if not
    int ccount;                    // size of copy message in bins
    int *cbinsrc;                  // list of my bin indices to copy from
    int *cbindest;                 // list of my bin indices to copy to
    MPI_Request *request;          // MPI datums for each recv message
    MPI_Status *status;
  };

  Swap *allswap;
  Swap *reverseswap;
  Swap **sectorswap;
  Swap **sectorreverseswap;
  int nsector;

  class AppOffLattice *appoff;

  int sitecustom;
  int ninteger,ndouble;

  double xprd,yprd,zprd;

  Swap *create_swap_all();
  Swap *create_swap_all_reverse();
  Swap *create_swap_sector(int, int *);
  Swap *create_swap_sector_reverse(int, int *);
  void free_swap(Swap *);

  void perform_swap(Swap *);
};

}
