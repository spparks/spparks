/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_NFW_H
#define APP_GRAIN_NFW_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrainNfw : public App {
 public:
  AppGrainNfw(class SPK *, int, char **);
  virtual ~AppGrainNfw();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);

 protected:
  int me,nprocs;
  int dimension;
  int nx_global,ny_global;            // size of global lattice [0,Nglobal-1]
                                      // global lattice has no ghosts
  int nx_local,ny_local;              // size of local lattice [0,Nlocal+1]
                                      // local lattice includes ghosts
  int nx_half,ny_half;                // indices of my local cell that is
                                      //   lower-left within upper-right quad
  int nx_offset,ny_offset;            // global indices (0:Nglobal-1) of
                                      //   my lower-left owned cell (1,1)
                                      // So, for single proc, offsets are zero
  int **lattice;                      // local lattice with ghost cells
                                      //   allocated to Nlocal+2 in each dim
                                      //   owned cells indexed 1 to Nlocal
                                      //   ghost cells = 0 and Nlocal+1

  int nx_procs,ny_procs;              // # of procs in each dim
  int procwest,proceast,procsouth,procnorth;    // neighbor procs

  int nspins;                          // # of possible spins
  int seed;                            // random number generator seed
  int ntimestep;
  double time,stoptime;
  FILE *fp;
  double temperature;
  int *buf;
  int maxbuf;
  double stats_time,stats_delta;       // statistics
  double dump_time,dump_delta;  

  double* propensity;

  class RandomPark *random;
  class CommGrain2D *comm;

  enum {nsector = 4};
  struct {
    int xlo,xhi,ylo,yhi;      // inclusive start/stop indices in each quadrant
  } quad[nsector];

  void iterate();
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

  void update_spin(const int&, const int&);
  void update_propensity(const int&, const int&) const;
  void init_propensity();
  double compute_propensity(const int&, const int&) const;
  void survey_neighbor(const int&, const int&, int&, int[], int[]) const;
};

}

#endif
