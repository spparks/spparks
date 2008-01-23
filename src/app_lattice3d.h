/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_LATTICE3D_H
#define APP_LATTICE3D_H

#include "stdio.h"
#include "app.h"

#include <stack>
#include "cluster.h"

namespace SPPARKS {

class AppLattice3d : public App {
  friend class SweepLattice3d;
  friend class DiagCluster3d;
  friend class DiagEprof3d;
  friend class DiagEnergy3d;

 public:
  AppLattice3d(class SPK *, int, char **);
  virtual ~AppLattice3d();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  virtual double site_energy(int, int, int) = 0;
  virtual int site_pick_random(int, int, int, double) = 0;
  virtual int site_pick_local(int, int, int, double) = 0;
  virtual double site_propensity(int, int, int, int) = 0;
  virtual void site_event(int, int, int, int) = 0;
  virtual void site_update_ghosts(int, int, int) = 0;
  virtual void site_clear_mask(char ***, int, int, int) = 0;

 protected:
  enum InitStyles {RANDOM,READ};
  int me,nprocs;
  int ntimestep,seed;
  int dump_style;
  int init_style;
  double time,stoptime;
  double temperature,t_inverse;
  int nsweep;

  int nx_global,ny_global,nz_global;     // global lattice (0 to nglobal-1)
  int nx_local,ny_local,nz_local ;       // local lattice (1 to nlocal)
                                         // does not include ghost sites
  int nx_offset,ny_offset,nz_offset;     // global indices of my (1,1) site
  int nx_sector_lo,nx_sector_hi;         // bounds of current sector
  int ny_sector_lo,ny_sector_hi;         // as set by sweeper
  int nz_sector_lo,nz_sector_hi;
  int nyz_local;

  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi; // upper/lower limits for local lattice
                                     // w/ ghost layer of thickness = delghost
                                     // local sites are from 1 to nlocal

  int ***lattice;                 // owned sites + ghost sites
  double *propensity;             // probability for each owned site
  int ***ijk2site;                // mapping of owned lattice to site index
  int **site2ijk;                 // mapping of owned sites to lattice indices

  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition
  int procwest,proceast;            // my neighbor procs
  int procsouth,procnorth;
  int procdown,procup;

  double masklimit;            // app-specific, used by sweeper

  int delghost,dellocal;       // app-specific thickness of 
                               // ghost and local layers needed for comm
                               // delghost affects upper and lower
                               // limits for local lattice
  FILE *fp;
  int *ibufdump, *ibufread;
  double *dbufdump;
  int maxdumpbuf;

  class RandomPark *random;
  class CommLattice3d *comm;

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}

  void iterate();
  void stats(char *);
  void stats_header(char *);
  void dump_header();
  void dump();
  void dump_lattice();
  void dump_coord();
  virtual void box_bounds(double *, double *, double *,
			  double *, double *, double *);
  virtual void xyz(int, int, int, double *, double *, double *);

  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void procs2lattice();

  void ijkpbc(int &, int &, int &);

  void read_spins(const char*);

  virtual void push_connected_neighbors(int, int, int, int***, int, std::stack<int>*);
  virtual void connected_ghosts(int, int, int, int***, Cluster*, int);

};

// remap i,j,k indices via PBC if needed

inline void AppLattice3d::ijkpbc(int &i, int &j, int &k)
{
  if (i < 1) i += nx_local;
  else if (i > nx_local) i -= nx_local;
  if (j < 1) j += ny_local;
  else if (j > ny_local) j -= ny_local;
  if (k < 1) k += nz_local;
  else if (k > nz_local) k -= nz_local;
}

}

#endif
