/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_LATTICE3D_H
#define APP_LATTICE3D_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppLattice3d : public App {
 public:
  AppLattice3d(class SPK *, int, char **);
  virtual ~AppLattice3d();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  double virtual site_energy(int, int, int) = 0;
  int virtual site_pick_random(int, int, int, double) = 0;
  int virtual site_pick_local(int, int, int, double) = 0;
  double virtual site_propensity(int, int, int) = 0;
  void virtual site_event(int, int, int) = 0;
  void virtual site_update_ghost(int, int, int) = 0;
  void virtual site_clear_mask(char ***, int, int, int) = 0;

 protected:
  int me,nprocs;
  int ntimestep,seed;
  double time,stoptime;
  double stats_time,stats_delta;
  double dump_time,dump_delta;
  double temperature,t_inverse;
  int nsweep;

  int nx_global,ny_global,nz_global;     // global lattice (0 to nglobal-1)
  int nx_local,ny_local,nz_local ;       // local lattice (1 to nlocal)
                                         // 0 and nlocal+1 are ghosts
  int nx_offset,ny_offset,nz_offset;     // global indices of my (1,1) site
  int nyz_local;

  int ***lattice;                   // owned sites + ghost sites
  double *propensity;               // probability for each owned site

  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition
  int procwest,proceast;            // my neighbor procs
  int procsouth,procnorth;
  int procdown,procup;

  FILE *fp;
  int *dumpbuf;
  int maxdumpbuf;

  class RandomPark *random;
  class CommLattice3d *comm;

  void virtual input_app(char *, int, char **);
  void iterate();
  void stats();
  void dump_header();
  void dump();

  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void procs2lattice();

  void site2ijk(const int, int &, int &, int &) const;
  int ijk2site(const int, const int, const int) const;
  void ijkpbc(const int, const int, const int, int &, int &, int &);
};

// convert site (0 to N-1) to i,j,k indices (1 to N)

inline void AppLattice3d::site2ijk(const int isite,
				   int &i, int &j, int &k) const
{
  i = isite/nyz_local + 1;
  j = isite/nz_local % ny_local + 1;
  k = isite % nz_local + 1;
}

// convert i,j,k indices (1 to N) to site (0 to N-1)

inline int AppLattice3d::ijk2site(const int i, const int j, const int k) const
{
  return (i-1)*nyz_local+(j-1)*nz_local+k-1;
}

// remap i,j,k indices (0 to N+1) to ii,jj,kk indices (1 to N) via PBC

inline void AppLattice3d::ijkpbc(const int i, const int j, const int k,
				 int &ii, int &jj, int &kk)
{
  if (i < 1) ii = i + nx_local;
  else if (i > nx_local) ii = i - nx_local;
  else ii = i;
  if (j < 1) jj = j + ny_local;
  else if (j > ny_local) jj = j - ny_local;
  else jj = j;
  if (k < 1) kk = k + nz_local;
  else if (k > nz_local) kk = k - nz_local;
  else kk = k;
}

}

#endif
