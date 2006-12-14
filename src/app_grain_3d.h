/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_3D_H
#define APP_GRAIN_3D_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrain3D : public App {
 public:
  AppGrain3D(class SPK *, int, char **);
  virtual ~AppGrain3D();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);

 protected:
  int me,nprocs;
  int ntimestep;
  int dimension;
  int nx_global,ny_global,nz_global;  // size of global lattice [0,Nglobal-1]
                                      // global lattice has no ghosts
  int nx_local,ny_local,nz_local;     // size of local lattice [0,Nlocal+1]
                                      // local lattice includes ghosts
  int nx_half,ny_half,nz_half;        // indices of my local cell that is
                                      //   lower-left within upper-right quad
  int nx_offset,ny_offset,nz_offset;  // global indices (0:Nglobal-1) of
                                      //   my lower-left owned cell (1,1)
                                      // So, for single proc, offsets are zero
  int ***lattice;                     // local lattice with ghost cells
                                      //   allocated to Nlocal+2 in each dim
                                      //   owned cells indexed 1 to Nlocal
                                      //   ghost cells = 0 and Nlocal+1

  int nx_procs,ny_procs,nz_procs;     // # of procs in each dim
  int procwest,proceast,procsouth,procnorth;    // neighbor procs
  int procdown,procup;    

  int nspins;                          // # of possible spins
  int nsweep;                          // # of sweeps to perform
  int seed;                            // random number generator seed
  int nstats,stats_next;
  int ndump,dump_next;
  FILE *fp;
  double temperature;
  int *buf;
  int maxbuf;

  class RandomPark *random;

  class CommGrain3D *comm;

  enum {nsector = 8};
  struct {
    int xlo,xhi,ylo,yhi,zlo,zhi;      // inclusive start/stop indices in each quadrant
  } quad[nsector];

  void iterate();
  void sweep(int);
  void stats();
  void dump_header();
  void dump();
  void dump_detailed(char*);

  void procs2lattice();
  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);
  double compute_energy();
  double energy_quadrant(int);
};

}

#endif
