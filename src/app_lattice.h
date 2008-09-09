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

#ifndef APP_LATTICE_H
#define APP_LATTICE_H

#include "stdio.h"
#include "app.h"

#include <stack>
#include "cluster.h"

namespace SPPARKS_NS {

class AppLattice : public App {
  friend class SweepLattice;
  friend class CommLattice;
  friend class DiagEnergy;
  friend class DiagCluster;

 public:
  AppLattice(class SPPARKS *, int, char **);
  virtual ~AppLattice();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  // pure virtual functions, must be defined in child class

  virtual double site_energy(int) = 0;
  virtual void site_event_rejection(int, class RandomPark *) = 0;
  virtual double site_propensity(int) = 0;
  virtual void site_event(int, class RandomPark *) = 0;

  void push_connected_neighbors(int, int*, int, std::stack<int>*);
  void connected_ghosts(int, int*, Cluster*, int);

  enum{NONE,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       RANDOM_2D,RANDOM_3D,FILENAME};

 protected:
  int me,nprocs;
  int ntimestep;
  double time,stoptime;
  double temperature,t_inverse;
  int nsweep;

  int latstyle;
  double latconst;
  int dimension;
  int nx,ny,nz;
  int nrandom;
  int latseed;
  double cutoff;
  char *latfile;

  int delpropensity;           // # of sites away needed to compute propensity
  int delevent;                // # of sites away affected by an event
  int numrandom;               // # of RN used by rejection routine

  bool Lmask;                  // from sweeper
  char *mask;

  double xprd,yprd,zprd;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // simulation box bounds
  double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;    // my portion of box

  int nglobal;                 // global # of sites
  int nlocal;                  // # of sites I own
  int nghost;                  // # of ghost sites I store

                               // these arrays stored for owned + ghost sites
  int *id;                     // global ID (1-N) of site
  int *owner;                  // proc who owns the site
  int *index;                  // index of site on owning proc
  double **xyz;                // coords of site

                               // per-site storage for owned + ghost sites
  int sitecustom;              // 0/1 for default or customized
  int ninteger,ndouble;        // # of int/double per site, 0,0 = just lattice
  int *lattice;                // default = single int value
  int **iarray;                // one or more ints per site
  double **darray;             // one or more doubles per site

  double *propensity;          // probabilities for each owned site
  int *i2site;                 // mapping of owned lattice to site index
  int *site2i;                 // mapping of owned sites to lattice index

                               // neigh info for owned sites
                               // and ghost sites up to delpropensity-1 layers
  int maxneigh;                // max neighbors of any site in entire system
  int *numneigh;               // # of neighbors of each site
  int **neighbor;              // list of neighbors of each site
                               // neighbor[i][j] =
                               // local index of Jth neigh of Ith owned site
                               // can point to owned or ghost site

  int nbasis;                  // basis atoms for regular lattices
  int ***cmap;                 // connectivity map for regular lattices

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition

  FILE *fp;
  FILE *fpdump;
  double *dbuf;
  int maxdumpbuf;

  enum DumpStyles {COORD,OPENDX};
  int dump_style;
  char* opendxroot;
  int opendxcount;

  class CommLattice *comm;
  class RandomPark *random;

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}

  void options(int, char **);
  void create_lattice();

  void structured_lattice();
  void random_lattice();
  void file_lattice();

  void ghosts_from_connectivity();
  void connectivity_within_cutoff();

  void iterate();
  void stats(char *);
  void stats_header(char *);
  void dump_header();
  void dump();

  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);

  void procs2lattice_2d();
  void procs2lattice_3d();

  int connect(int, int);
  void offsets();
};

}

#endif
