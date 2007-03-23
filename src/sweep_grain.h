/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_GRAIN_H
#define SWEEP_GRAIN_H

#include "sweep.h"

namespace SPPARKS {

class SweepGrain : public Sweep {
 public:
  SweepGrain(class SPK *, int, char **);
  ~SweepGrain();
  void input(int, char **){};
  void init(AppGrain*, const int, const int, const int, 
		      const int, const int, const int, 
		      const int, const int, const int, 
		      const int, const int, const int, 
		      int***, 
		      const int, const int, 
	              const int, const int, 
	              const int, const int, 
	    const int, const double, 
	    int (AppGrain::*)(int,int,int,int,int),
	    void (AppGrain::*)(char***,int,int,int,int) );
  void do_sweep();
  double compute_energy();
  double energy_quadrant(const int);

 private:
  void sweep_quadrant_2d(int);
  void sweep_quadrant_mask_2d(int);
  void sweep_quadrant_strict_2d(int, int);
  void sweep_quadrant_mask_strict_2d(int, int);

  void sweep_quadrant_3d(int);
  void sweep_quadrant_mask_3d(int);
  void sweep_quadrant_strict_3d(int, int);
  void sweep_quadrant_mask_strict_3d(int, int);

 private:
  AppGrain* appgrain;
  int me,nprocs;
  int dimension;
  int nx_global,ny_global,nz_global;  // size of global lattice [0,Nglobal-1]
                                      // global lattice has no ghosts
  int nx_local,ny_local,nz_local;     // size of local lattice [0,Nlocal+1]
                                      // local lattice includes ghosts
  int nx_offset,ny_offset,nz_offset;  // global indices (0:Nglobal-1) of
                                      //   my lower-left owned cell (1,1)
                                      // So, for single proc, offsets are zero
  int nx_half,ny_half,nz_half;        // indices of my local cell that is
                                      //   lower-left within upper-right quad
  int ***lattice;                     // local lattice with ghost cells
                                      //   allocated to Nlocal+2 in each dim
                                      //   owned cells indexed 1 to Nlocal
                                      //   ghost cells = 0 and Nlocal+1
  int **lat_2d, ***lat_3d;            // 2D and 3D pointers to lattice

  int procwest,proceast;              // neighbor procs
  int procsouth,procnorth;
  int procdown,procup;    

  int nspins;                         // # of possible spins
  int seed;                           // random number generator seed
  double temperature;

  class RandomPark *random;
  class CommGrain2D *comm_2d;         // Pointer to 2D comm object
  class CommGrain3D *comm_3d;         // Pointer to 3D comm object

  // Data for quadrants
  int nquad;
  enum {nquadmax = 8};
  struct {
    int xlo,xhi,ylo,yhi,zhi,zlo;      // inclusive start/stop indices in each quadrant
  } quad[nquadmax];

  // Data for masking methods
  bool Lmask;                         // Flag indicating whether masking used
  int masklimit;                      // Limit value for for not masking
  char*** mask;                        // Local lattice of mask flags

  // Data for strict methods
  bool Lstrict;                       // Flag indicating whether sweeping is strictly
                                      // independent of number of processors.
  class RandomPark ***ranlat;          // array of random number generator pointers
  int ncolor;                         // Number of colors to use

  // Pointer-to-member functions for different lattice stencils
  int (AppGrain::*es_fp)(int,int,int,int,int);
  void (AppGrain::*um_fp)(char***,int,int,int,int);

};

}

#endif
