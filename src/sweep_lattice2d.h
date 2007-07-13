/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE2D_H
#define SWEEP_LATTICE2D_H

#include "sweep.h"

namespace SPPARKS {

class SweepLattice2d : public Sweep {
 public:
  SweepLattice2d(class SPK *, int, char **);
  ~SweepLattice2d();
  void init();
  void do_sweep(double &);
  double compute_energy() {return 0.0;};

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict;
  double delt;

  int nx_local,ny_local;
  int nx_offset,ny_offset;
  int **lattice;
  double temperature,t_inverse;

  class AppLattice2d *applattice;       
  class CommLattice2d *comm;
  class RandomPark *random;

  int ncolor;
  double masklimit;
  char **mask;
  class RandomPark **ranlat;

  int nquad;
  struct {
    int xlo,xhi,ylo,yhi;     // inclusive start/stop indices in each quadrant
  } quad[4];

  typedef void (SweepLattice2d::*FnPtr)(int, int);  // pointer to sweep method
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
