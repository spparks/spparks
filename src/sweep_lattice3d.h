/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE3D_H
#define SWEEP_LATTICE3D_H

#include "sweep.h"

namespace SPPARKS {

class SweepLattice3d : public Sweep {
 public:
  SweepLattice3d(class SPK *, int, char **);
  ~SweepLattice3d();
  void init(const int, const int, const int,
	    const int, const int, const int,
	    const int, const int, const int,
	    const int, const int, const int, const int, const int, const int,
	    int ***, const double);
  void do_sweep(double &);
  double compute_energy() {return 0.0;};

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict;
  double delt;

  int nx_local,ny_local,nz_local;
  int nx_offset,ny_offset,nz_offset;
  int ***lattice;
  double temperature,t_inverse;

  class AppLattice3d *applattice;       
  class CommLattice3d *comm;
  class RandomPark *random;

  int ncolor;
  double masklimit;
  char ***mask;
  class RandomPark ***ranlat;

  int nquad;
  struct {
    int xlo,xhi,ylo,yhi,zlo,zhi; // inclusive start/stop indices in each octant
  } quad[8];

  typedef void (SweepLattice3d::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sector;                   

  // sweep methods

  void sweep_quadrant(int, int);
  void sweep_quadrant_mask(int, int);
  void sweep_quadrant_mask_picklocal(int, int);
  void sweep_quadrant_strict(int, int);
  void sweep_quadrant_mask_strict(int, int);
};

}

#endif
