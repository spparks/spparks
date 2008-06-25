/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE2D_H
#define SWEEP_LATTICE2D_H

#include "sweep.h"

namespace SPPARKS_NS {

class SweepLattice2d : public Sweep {
  friend class AppLattice2d;

 public:
  SweepLattice2d(class SPPARKS *, int, char **);
  ~SweepLattice2d();
  void init();
  void do_sweep(double &);

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict,Lkmc,Ladapt;
  double delt;

  int nx_local,ny_local;
  int nx_offset,ny_offset;
  int **lattice,**ij2site;
  double temperature,t_inverse;

  class AppLattice2d *applattice;       
  class CommLattice2d *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char **mask;
  class RandomPark **ranlat;

  double masklimit;              // App-specific settings
  int delghost,dellocal;
  int nxlo,nxhi,nylo,nyhi;

  int nsector;
  struct {
    int xlo,xhi,ylo,yhi;     // inclusive start/stop indices in this quadrant
    int nx,ny;               // size of quadrant
    class Solve *solve;      // KMC solver
    double *propensity;      // propensities for quadrant sites
    int **site2ij;           // map from quadrant sites to lattice index
    int *sites;              // list of sites to pass to solver
  } sector[4];

  typedef void (SweepLattice2d::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  double pmax,pmaxall,deln0;
  // sweep methods

  void sweep_sector(int, int);
  void sweep_sector_mask(int, int);
  void sweep_sector_mask_picklocal(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  void boundary_clear_mask(int, int, int, int);
};

}

#endif
