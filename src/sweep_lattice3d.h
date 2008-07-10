/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE3D_H
#define SWEEP_LATTICE3D_H

#include "sweep.h"

namespace SPPARKS_NS {

class SweepLattice3d : public Sweep {
  friend class AppLattice3d;

 public:
  SweepLattice3d(class SPPARKS *, int, char **);
  ~SweepLattice3d();
  void init();
  void do_sweep(double &);

 private:
  int seed;
  double delt;

  int nx_local,ny_local,nz_local;
  int nx_offset,ny_offset,nz_offset;
  int ***lattice;
  int ***ijk2site;

  class AppLattice3d *applattice;       
  class CommLattice3d *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char ***mask;
  class RandomPark ***ranlat;

  int delpropensity,delevent;            // app-specific settings
  int numrandom;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;

  double pmax,pmaxall,deln0;

  int nsector;
  struct {
    int xlo,xhi,ylo,yhi,zlo,zhi; // inclusive start/stop indices in this octant
    int nx,ny,nz;                // size of octant
    int nborder;                 // # of sites with non-sector site as neighbor
    class Solve *solve;          // KMC solver
    int **site2ijk;              // map from octant sites to lattice index
    int ***ijk2site;             // map from lattice index to sector sites
    int **border;                // lattice index for each border site
    int *sites;                  // list of sites to pass to solver
    double *propensity;          // propensities for sector sites
  } sector[8];

  // sweep methods

  typedef void (SweepLattice3d::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  void sweep_sector(int, int);
  void sweep_sector_mask(int, int);
  void sweep_sector_mask_picklocal(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  int find_border_sites(int);
  void boundary_clear_mask(int);
};

}

#endif
