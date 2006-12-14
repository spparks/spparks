/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_GRAIN_3D_H
#define COMM_GRAIN_3D_H

#include "sysptr.h"

namespace SPPARKS {

class CommGrain3D : protected SysPtr  {
 public:
  explicit CommGrain3D(class SPK *);
  ~CommGrain3D();
  void setup(const int, const int, const int, const int, 
             const int, const int, const int, const int, 
	     const int, const int, const int, const int);
  void communicate(int ***, const int);

 private:
  int *sendbuf;          // send buffer for all comm
  int *recvbuf;          // recv buffer for all comm
  int maxsend;           // current size of send buffer
  int maxrecv;           // current size of recv buffer

  int nswap;             // # of swaps to perform, set to 2 for 2D Potts model
  int maxswap;           // max # of swaps memory is allocated for
  int nsector;           // # of sectors, set to 4 for 2D Potts model
  int maxsector;         // max # of sectors memory is allocated for

  struct SwapInfo {
    int numx, numy, numz;    // # of data to send/recv in each swap
    int sendproc, recvproc;  // proc to send/recv to/from at each swap
    int sendix, recvix;      // Location of first data value for each send/re cv
    int sendiy, recviy;      // ix, iy and iz indexes for 3D data array
    int sendiz, recviz;      // 

    int copyixa, copyixb;    // Copy from/to locations for fragments
    int copyiya, copyiyb;    // a,b refer to locations of original, copy
                             // swap 0: iz comes from recviz/sendiz
                             // swap 1: only ix needed

  };

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects

  void allocate_swap(const int, const int);          // allocate swap arrays
  void free_swap();                                  // free swap arrays
};

}

#endif
