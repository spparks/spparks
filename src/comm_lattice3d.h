/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_LATTICE3D_H
#define COMM_LATTICE3D_H

#include "sysptr.h"

namespace SPPARKS {

class CommLattice3d : protected SysPtr {
 public:
  CommLattice3d(class SPK *);
  ~CommLattice3d();
  void init(const int, const int, const int,
	    const int, const int, const int,
	    const int, const int, const int);
  void sector(int ***, const int);
  void all(int ***);

 private:
  int me,nprocs;
  int nx_local,ny_local,nz_local;
  int procwest,proceast,procsouth,procnorth,procdown,procup;

  int *sendbuf;          // send buffer for all comm
  int *recvbuf;          // recv buffer for all comm
  int maxbuf;            // size of buffers

  int nswap;             // # of swaps to perform per sector
  int maxswap;           // max # of swaps memory is allocated for
  int nsector;           // # of sectors
  int maxsector;         // max # of sectors memory is allocated for

  struct SwapInfo {
    int numx, numy, numz;    // # of data to send/recv in each swap = numx*numy*numz
    int sendproc, recvproc;  // proc to send/recv to/from at each swap
    int sendix, recvix;      // loc of first data value for each send/recv
    int sendiy, recviy;      // ix, iy and iz indexes for 3D data array
    int sendiz, recviz;
    int copyixa, copyixb;    // copy from/to locations for fragments
    int copyiya, copyiyb;    // a,b refer to locations of original
                             // swap 0: iz comes from recviz/sendiz
                             // swap 1: only ix needed
                             // swap 2: no copy needed
  };

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects

  void allocate_swap(const int, const int);          // allocate swap arrays
  void free_swap();                                  // free swap arrays
};

}

#endif

