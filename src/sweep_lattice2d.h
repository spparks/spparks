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
  double delt;

  int nx_local,ny_local;
  int nx_offset,ny_offset;
  int **lattice;
  int **ij2site;

  class AppLattice2d *applattice;       
  class CommLattice2d *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char **mask;
  class RandomPark **ranlat;

  int delpropensity,delevent;     // app-specific settings
  int numrandom;
  int nxlo,nxhi,nylo,nyhi;

  double pmax,pmaxall,deln0;

  int nsector;
  struct {
    int xlo,xhi,ylo,yhi;     // inclusive start/stop indices in this quadrant
    int nx,ny;               // size of quadrant
    int nborder;             // # of sites with non-sector site as neighbor
    class Solve *solve;      // KMC solver
    int **site2ij;           // map from sector sites to lattice index
    int **ij2site;           // map from lattice index to sector sites
    int **border;            // lattice index for each border site
    int *sites;              // list of sites to pass to solver
    double *propensity;      // propensities for sector sites
  } sector[4];

  // sweep methods

  typedef void (SweepLattice2d::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  void sweep_sector(int, int);
  void sweep_sector_mask(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  int find_border_sites(int);
  void boundary_clear_mask(int);
};

}

#endif
