/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_LATTICE2D_H
#define APP_LATTICE2D_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppLattice2d : public App {
  friend class SweepLattice2d;

 public:
  AppLattice2d(class SPK *, int, char **);
  virtual ~AppLattice2d();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  double virtual site_energy(int, int) = 0;
  int virtual site_pick_random(int, int, double) = 0;
  int virtual site_pick_local(int, int, double) = 0;
  double virtual site_propensity(int, int, int) = 0;
  void virtual site_event(int, int, int) = 0;
  void virtual site_update_ghosts(int, int) = 0;
  void virtual site_clear_mask(char **, int, int) = 0;

 protected:
  int me,nprocs;
  int ntimestep,seed;
  double time,stoptime;
  double stats_time,stats_delta;
  double dump_time,dump_delta;
  double temperature,t_inverse;
  int nsweep;

  int nx_global,ny_global;           // global lattice (0 to nglobal-1)
  int nx_local,ny_local;             // local lattice (1 to nlocal)
                                     // 0 and nlocal+1 are ghosts
  int nx_offset,ny_offset;           // global indices of my (1,1) site
  int nx_sector_lo,nx_sector_hi;     // bounds of current sector
  int ny_sector_lo,ny_sector_hi;     // as set by sweeper

  int **lattice;               // owned lattice + ghost lattice
  double *propensity;          // probability for each owned site
  int **ij2site;               // mapping of owned lattice to sites
  int **site2ij;               // mapping of owned sites to lattice indices

  int nx_procs,ny_procs;       // procs in each dim of lattice partition
  int procwest,proceast;       // my neighbor procs
  int procsouth,procnorth;

  double masklimit;            // app-specific, used by sweeper

  FILE *fp;
  int *dumpbuf;
  int maxdumpbuf;

  class RandomPark *random;
  class CommLattice2d *comm;

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}

  void iterate();
  void stats();
  void dump_header();
  void dump();

  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void procs2lattice();
  void ijpbc(int &, int &);
};

// remap i,j indices via PBC if needed

inline void AppLattice2d::ijpbc(int &i, int &j)
{
  if (i < 1) i += nx_local;
  else if (i > nx_local) i -= nx_local;
  if (j < 1) j += ny_local;
  else if (j > ny_local) j -= ny_local;
}

}

#endif
