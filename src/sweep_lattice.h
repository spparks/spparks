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
  friend class CommLattice;

 public:
  SweepLattice(class SPK *, int, char **);
  ~SweepLattice();
  void init();
  void do_sweep(double &);

 private:
  int seed;
  bool Lmask,Lpicklocal,Lstrict,Lkmc,Ladapt;
  double delt;

  int *lattice;

  int dimension;
  int nlocal;
  int *i2site;
  double temperature,t_inverse;

  class AppLattice *applattice;       
  class CommLattice *comm;
  class RandomPark *random;

  int ncolor,delcol;
  char *mask;
  class RandomPark *ranlat;

  double masklimit;              // App-specific settings
  int delghost,dellocal;

  int nsector;
  struct {
    int nlocal;              // # of owned sites in sector
    int nmax;                // max # of sites vectors can hold
    int nborder;             // # of owned sites with a ghost site as neighbor
    class Solve *solve;      // KMC solver
    double *propensity;      // propensities for sector sites
    int *i2site;             // map from lattice index to sector sites
    int *site2i;             // map from sector sites to lattice index
    int *sites;              // list of sites to pass to solver
    int *border;             // lattice index for each border site
  } sector[8];

  int find_border_sites(int, int *, int *, int **, int **);

  typedef void (SweepLattice::*FnPtr)(int, int);  // pointer to sweep method
  FnPtr sweeper;

  double pmax,pmaxall,deln0;

  // sweep methods

  void sweep_sector_lattice(int, int);
  void sweep_sector_data(int, int);

  void sweep_sector_mask(int, int);
  void sweep_sector_mask_picklocal(int, int);
  void sweep_sector_strict(int, int);
  void sweep_sector_mask_strict(int, int);
  void sweep_sector_kmc(int, int);

  void boundary_clear_mask();
};

}

#endif
