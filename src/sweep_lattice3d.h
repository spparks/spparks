/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE3D_H
#define SWEEP_LATTICE3D_H

#include "sweep.h"

namespace SPPARKS {

class SweepLattice3d : public Sweep {
  friend class AppLattice3d;

 public:
  SweepLattice3d(class SPK *, int, char **);
  ~SweepLattice3d();
  void init();
  void do_sweep(double &);
  void boundary_clear_mask(int, int, int, int, int, int);

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict,Lkmc;
  double delt;

  int nx_local,ny_local,nz_local;
  int nx_offset,ny_offset,nz_offset;
  int ***lattice,***ijk2site;
  double temperature,t_inverse;

  class AppLattice3d *applattice;       
  class CommLattice3d *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char ***mask;
  class RandomPark ***ranlat;

  double masklimit;              // App-specific settings
  int delghost,dellocal;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;

  int nquad;
  struct {
    int xlo,xhi,ylo,yhi,zlo,zhi; // inclusive start/stop indices in this octant
    int nx,ny,nz;                // size of octant
    int nsites;                  // total sites in octant
    class Solve *solve;          // KMC solver
    double *propensity;          // propensities for octant sites
    int **site2ijk;              // map from octant sites to local lattice
    int *sites;                  // list of sites to pass to solver
  } quad[8];

  typedef void (SweepLattice3d::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sector;                   

  // sweep methods

  void sweep_quadrant(int, int);
  void sweep_quadrant_mask(int, int);
  void sweep_quadrant_mask_picklocal(int, int);
  void sweep_quadrant_strict(int, int);
  void sweep_quadrant_mask_strict(int, int);
  void sweep_quadrant_kmc(int, int);
};

}

#endif
