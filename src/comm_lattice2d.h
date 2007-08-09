/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_LATTICE2D_H
#define COMM_LATTICE2D_H

#include "sysptr.h"

namespace SPPARKS {

class CommLattice2d : protected SysPtr {
 public:
  CommLattice2d(class SPK *);
  ~CommLattice2d();
  void init(const int, const int, 
	    const int, const int, const int, const int,
	    const int, const int);
  void sector(int **, const int);
  void all(int **);

 private:
  int me,nprocs;
  int nx_local,ny_local;
  int procwest,proceast,procsouth,procnorth;
  int delghost,dellocal; // thickness of ghost and local communication layers 

  int *sendbuf;          // send buffer for all comm
  int *recvbuf;          // recv buffer for all comm
  int maxbuf;            // size of buffers

  int nswap;             // # of swaps to perform per sector
  int maxswap;           // max # of swaps memory is allocated for
  int nsector;           // # of sectors
  int maxsector;         // max # of sectors memory is allocated for

  struct SwapInfo {
    int stride;              // spacing of data values in memory for this swap
    int numx, numy;          // # of data to send/recv in each swap = numx*numy
    int sendproc, recvproc;  // proc to send/recv to/from at each swap
    int sendix, recvix;      // loc of first data value for each send/recv
    int sendiy, recviy;      // ix and iy indexes for 2D data array
    int copyixa, copyixb;    // copy from/to locations for corners
                             // a,b refer to locations of original
                             // swap 0: iy comes from recviy/sendiy
                             // swap 1: no copy needed
  };

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects

  void allocate_swap(const int, const int);          // allocate swap arrays
  void free_swap();                                  // free swap arrays
};

}

#endif

