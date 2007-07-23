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
  double virtual site_propensity(int, int) = 0;
  void virtual site_event(int, int) = 0;
  // This is a short-cut; need to make this pure virtual
  void virtual site_event_sector(int, int) {};
  void virtual site_update_ghost(int, int) = 0;
  void virtual site_clear_mask(char **, int, int) = 0;

 protected:
  int me,nprocs;
  int ntimestep,seed;
  double time,stoptime;
  double stats_time,stats_delta;
  double dump_time,dump_delta;
  double temperature,t_inverse;
  int nsweep;

  int nx_global,ny_global;     // global lattice (0 to nglobal-1)
  int nx_local,ny_local;       // local lattice (1 to nlocal)
                               // 0 and nlocal+1 are ghosts
  int nx_offset,ny_offset;     // global indices of my (1,1) site

  int **lattice;               // owned sites + ghost sites
  double *propensity;          // probability for each owned site

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

  void site2ij(const int, int &, int &) const;
  int ij2site(const int, const int) const;
  void ijpbc(const int, const int, int &, int &);
};

// convert site (0 to N-1) to i,j indices (1 to N)

inline void AppLattice2d::site2ij(const int isite, int &i, int &j) const
{
  i = (isite / ny_local) + 1;
  j = isite - (i-1)*ny_local + 1;
}

// convert i,j indices (1 to N) to site (0 to N-1)

inline int AppLattice2d::ij2site(const int i, const int j) const
{
  return (i-1)*ny_local + j-1;
}

// remap i,j indices (0 to N+1) to ii,jj indices (1 to N) via PBC

inline void AppLattice2d::ijpbc(const int i, const int j, int &ii, int &jj)
{
  if (i < 1) ii = i + nx_local;
  else if (i > nx_local) ii = i - nx_local;
  else ii = i;
  if (j < 1) jj = j + ny_local;
  else if (j > ny_local) jj = j - ny_local;
  else jj = j;
}

}

#endif
