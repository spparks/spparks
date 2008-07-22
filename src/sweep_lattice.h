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

#ifndef SWEEP_LATTICE_H
#define SWEEP_LATTICE_H

#include "sweep.h"

namespace SPPARKS_NS {

class SweepLattice : public Sweep {
  friend class AppLattice;
  friend class CommLattice;

 public:
  SweepLattice(class SPPARKS *, int, char **);
  ~SweepLattice();
  void init();
  void do_sweep(double &);

 private:
  int seed;
  double delt;

  int dimension;
  int nlocal;
  int *i2site;

  class AppLattice *applattice;       
  class CommLattice *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char *mask;
  class RandomPark *ranlat;

  double pmax,pmaxall,deln0;

  int nsector;
  struct {
    int nlocal;              // # of owned sites in sector
    int nmax;                // max # of sites vectors can hold
    int nborder;             // # of sites with non-sector site as neighbor
    class Solve *solve;      // KMC solver
    int *site2i;             // map from sector sites to lattice index
    int *i2site;             // map from lattice index to sector sites
    int *border;             // lattice index for each border site
    int *sites;              // list of sites to pass to solver
    double *propensity;      // propensities for sector sites
  } sector[8];

  // sweep methods

  typedef void (SweepLattice::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  void sweep_sector(int, int);
  void sweep_sector_mask(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  int find_border_sites(int, int, int, int *, int **);
  void boundary_clear_mask(int);
};

}

#endif
