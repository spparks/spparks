/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_H
#define APP_GRAIN_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrain : public App {
 public:
  AppGrain(class SPK *, int, char **);
  virtual ~AppGrain();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);
  // Energy functions for different lattice stencils
  // These must all look the same:
  // int energy energy_site_square_8nn(const int ik, const int i, const int j, const int k, const int ib);
  // energy = site energy
  // i,j,k are the site indexes
  // ib is a basis index (could be used for things like triangular lattice) 
  int energy_site_square_8nn(const int, const int, const int, const int, const int);
  int energy_site_square_4nn(const int, const int, const int, const int, const int);
  int energy_site_cubic_26nn(const int, const int, const int, const int, const int);
  int energy_site_cubic_6nn(const int, const int, const int, const int, const int);
  // Mask functions for different lattice stencils
  void update_mask_square_8nn(char***, const int, const int, const int, const int);
  void update_mask_square_4nn(char***, const int, const int, const int, const int);
  void update_mask_cubic_26nn(char***, const int, const int, const int, const int);
  void update_mask_cubic_6nn(char***, const int, const int, const int, const int);
  int pick_local_square_8nn(const int, const int, const int, const int, const double);

  // N-Fold Way functions for different lattice stencils
  void flip_spin_square_8nn(const int&, const int&, const int&, const int&);
  void flip_spin_cubic_26nn(const int&, const int&, const int&, const int&);
  void update_propensity_square_8nn(const int&, const int&, const int&, const int&) const;
  void update_propensity_cubic_26nn(const int&, const int&, const int&, const int&) const;
  double compute_propensity_square_8nn(const int&, const int&, const int&, const int&) const;
  double compute_propensity_cubic_26nn(const int&, const int&, const int&, const int&) const;
  
 protected:
  int me,nprocs;
  int ntimestep;
  int dimension;
  double time,stoptime;
  double stats_time,stats_delta;       // statistics
  double dump_time,dump_delta;         // dump
  int sweep_create;                    // set if App creates its own Sweeper

  double* propensity;

  int nx_global,ny_global,nz_global;  // size of global lattice [0,Nglobal-1]
                                      // global lattice has no ghosts
  int nx_local,ny_local,nz_local;     // size of local lattice [0,Nlocal+1]
                                      // local lattice includes ghosts
  int nyz_local;                      // Frequently used product
  int nx_offset,ny_offset,nz_offset;  // global indices (0:Nglobal-1) of
                                      //   my lower-left owned cell (1,1)
                                      // So, for single proc, offsets are zero
  int ***lattice;                     // local lattice with ghost cells
                                      //   allocated to Nlocal+2 in each dim
                                      //   owned cells indexed 1 to Nlocal
                                      //   ghost cells = 0 and Nlocal+1
  int **lat_2d, ***lat_3d;            // 2D and 3D pointers to lattice

  int nx_procs,ny_procs,nz_procs;     // # of procs in each dim
  int procwest,proceast;               // neighbor procs
  int procsouth,procnorth;
  int procdown,procup;    
  int nspins;                          // # of possible spins
  int nsweep;                          // # of sweeps to perform
  int seed;                            // random number generator seed
  FILE *fp;
  double temperature;
  int *dumpbuf;
  int maxdumpbuf;
  // Pointer-to-member functions for different lattice stencils
  int (AppGrain::*es_fp)(int,int,int,int,int);
  void (AppGrain::*um_fp)(char***,int,int,int,int);
  int (AppGrain::*pl_fp)(int,int,int,int,double);
  void (AppGrain::*fs_fp)(const int&, const int&, const int&, const int&);
  double (AppGrain::*cp_fp)(const int&, const int&, const int&, const int&) const;
  void (AppGrain::*up_fp)(const int&, const int&, const int&, const int&) const;

  class RandomPark *random;

  void iterate();
  void stats();
  void dump_header();
  void dump();
  void dump_detailed(char*);

  void procs2lattice_2d();
  void procs2lattice_3d();
  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void init_propensity();
  void survey_neighbor(const int&, const int&, int&, int[], int[]) const;

  void isite2ij(const int&, int&, int&) const;
  int ij2isite(const int&, const int&) const;
  void isite2ijk(const int&, int&, int&, int&) const;
  int ijk2isite(const int&, const int&, const int&) const;

};

inline void AppGrain::isite2ij(const int& isite, int& i, int& j) const {
  i = (isite / ny_local) + 1;
  j = isite - (i-1)*ny_local + 1;
}

inline int AppGrain::ij2isite(const int& i, const int& j) const {
  return (i-1)*ny_local+j-1;
}

inline void AppGrain::isite2ijk(const int& isite, int& i, int& j, int& k) const {
  i = isite/nyz_local + 1;
  j = isite/nz_local % ny_local + 1;
  k = isite % nz_local + 1;
}

inline int AppGrain::ijk2isite(const int& i, const int& j, const int& k) const {
  return (i-1)*nyz_local+(j-1)*nz_local+k-1;
}


}

#endif
