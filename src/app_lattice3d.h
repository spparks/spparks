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

#ifndef APP_LATTICE3D_H
#define APP_LATTICE3D_H

#include "stdio.h"
#include "app.h"

#include <stack>
#include "cluster.h"

namespace SPPARKS_NS {

class AppLattice3d : public App {
  friend class SweepLattice3d;
  friend class DiagCluster3d;
  friend class DiagEprof3d;
  friend class DiagEnergy3d;

 public:
  AppLattice3d(class SPPARKS *, int, char **);
  virtual ~AppLattice3d();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  // pure virtual functions, must be defined in child class

  virtual double site_energy(int, int, int) = 0;
  virtual void site_event_rejection(int, int, int, class RandomPark *) = 0;
  virtual double site_propensity(int, int, int) = 0;
  virtual void site_event(int, int, int, int, class RandomPark *) = 0;

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

  int delpropensity;           // # of sites away needed to compute propensity
  int delevent;                // # of sites away affected by an event
  int numrandom;               // # of RN used by rejection routine

  bool Lmask;                  // from sweeper
  char ***mask;

  FILE *fp;
  int *ibufdump, *ibufread;
  double *dbufdump;
  int maxdumpbuf;

  class RandomPark *random;
  class CommLattice3d *comm;

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}

  void iterate();
  void update_ghost_sites(int, int, int);
  void add_unique(int, int &, int *, int, int, int);

  void stats(char *);
  void stats_header(char *);
  void dump_header();
  void dump();
  void dump_lattice();
  void dump_coord();
  void box_bounds(double *, double *, double *,
		  double *, double *, double *);
  void xyz(int, int, int, double *, double *, double *);

  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void procs2lattice();
  void ijkpbc(int &, int &, int &);

  void read_spins(const char*);

  virtual void push_connected_neighbors(int, int, int, int***, int,
					std::stack<int>*);
  virtual void connected_ghosts(int, int, int, int***, Cluster*, int);

};

// remap i,j,k indices via global PBC

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
