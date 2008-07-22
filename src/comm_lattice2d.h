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

#ifndef COMM_LATTICE2D_H
#define COMM_LATTICE2D_H

#include "pointers.h"

namespace SPPARKS_NS {

class CommLattice2d : protected Pointers {
 public: 
  CommLattice2d(class SPPARKS *);
  ~CommLattice2d();
  void init(const int, const int, 
	    const int, const int, const int, const int,
	    const int, const int);
  void all(int **);
  void sector(int **, const int);
  void reverse_sector(int **, const int);

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

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects (forward comm)
  SwapInfo** reverseinfo;    // nsector x nswap 2D array of SwapInfo objects (reverse comm)

  void setup_swapinfo();
  void setup_reverseinfo();
  void sector_multilayer(int **, const int);
  void sector_multilayer_destroy(int **, const int);
  void reverse_sector_multilayer(int **, const int);
  void all_multilayer(int **);
  void allocate_swap(const int, const int);          // allocate swap arrays
};

}

#endif

