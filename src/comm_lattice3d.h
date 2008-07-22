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

#ifndef COMM_LATTICE3D_H
#define COMM_LATTICE3D_H

#include "pointers.h"

namespace SPPARKS_NS {

class CommLattice3d : protected Pointers {
 public:
  CommLattice3d(class SPPARKS *);
  ~CommLattice3d();
  void init(const int, const int, const int,
	    const int, const int, const int,
	    const int, const int, const int,
	    const int, const int);
  void all(int ***);
  void sector(int ***, const int);
  void reverse_sector(int ***, const int);

 private:
  int me,nprocs;
  int nx_local,ny_local,nz_local;
  int procwest,proceast,procsouth,procnorth,procdown,procup;
  int delghost,dellocal; // thickness of ghost and local communication layers 

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

  SwapInfo** swapinfo;       // nsector x nswap 2D array of SwapInfo objects (forward comm)
  SwapInfo** reverseinfo;    // nsector x nswap 2D array of SwapInfo objects (reverse comm)

  void setup_swapinfo();
  void setup_reverseinfo();
  void sector_onelayer(int ***, const int);
  void sector_multilayer(int ***, const int);
  void sector_multilayer_destroy(int ***, const int);
  void all_onelayer(int ***);
  void all_multilayer(int ***);
  void reverse_sector_onelayer(int ***, const int);
  void reverse_sector_multilayer(int ***, const int);
  void reverse_all(int ***);
  void reverse_all_onelayer(int ***);
  void reverse_all_multilayer(int ***);
  void allocate_swap(const int, const int);          // allocate swap arrays
};

}

#endif
