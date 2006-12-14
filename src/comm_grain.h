/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_GRAIN_H
#define COMM_GRAIN_H

#include "sysptr.h"

namespace SPPARKS {

class CommGrain : protected SysPtr  {
 public:
  explicit CommGrain(class SPK *);
  ~CommGrain();
  void setup(const int, const int, const int, const int, 
	     const int, const int, const int, const int);
  void communicate(int **, const int);

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
    int stride;              // spacing of data values in memory for this swap
    int sendnum, recvnum;    // # of data to send/recv in each swap
    int sendproc, recvproc;  // proc to send/recv to/from at each swap
    int sendix, recvix;      // Location of first data value for each send/recv
    int sendiy, recviy;      // ix and iy indexes for 2D data array
    int copy0ixa, copy0ixb;  // From/To Locations for corner copy operation
                             // a,b refer to locations of original, copy
                             // ix indices; iy comes from recviy/sendiy)
  };

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects

  void allocate_swap(const int, const int);          // allocate swap arrays
  void free_swap();                                  // free swap arrays
};

}

#endif
