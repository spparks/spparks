/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWEEP_LATTICE_H
#define SWEEP_LATTICE_H

#include "sweep.h"

namespace SPPARKS {

class SweepLattice : public Sweep {
  friend class AppLattice;

 public:
  SweepLattice(class SPK *, int, char **);
  ~SweepLattice();
  void init();
  void do_sweep(double &);

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict,Lkmc;
  double delt;

  int *lattice;
  double temperature,t_inverse;

  class AppLattice *applattice;       
  class CommLattice *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char *mask;
  class RandomPark *ranlat;

  double masklimit;              // App-specific settings
  int delghost,dellocal;

  int dimension;
  int nlocal,nghost;

  int nsector;
  struct {
    int n;                   // size of sector
    class Solve *solve;      // KMC solver
    double *propensity;      // propensities for sector sites
  } sector[8];

  typedef void (SweepLattice::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  // sweep methods

  void sweep_sector(int, int);
  void sweep_sector_mask(int, int);
  void sweep_sector_mask_picklocal(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  void boundary_clear_mask();
};

}

#endif
