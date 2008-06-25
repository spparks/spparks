/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

  virtual double site_energy(int) = 0;
  virtual void site_pick_random(int ,double) = 0;
  virtual void site_pick_local(int, double) = 0;
  virtual double site_propensity(int, int) = 0;
  virtual void site_event(int, int) = 0;
  virtual void site_clear_mask(char *, int) = 0;

  void site_save(int);
  void site_restore(int);

  void push_connected_neighbors(int, int*, int, std::stack<int>*);
  void connected_ghosts(int, int*, Cluster*, int);

 protected:
  int me,nprocs;
  int ntimestep,seed;
  int dump_style;
  double time,stoptime;
  double temperature,t_inverse;
  int nsweep;

  int latstyle;
  double latconst;
  int dimension;
  int nx,ny,nz;
  int nrandom;
  double cutoff;
  char *latfile;

  int delghost,dellocal;
  int masklimit;

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

  int maxneigh;                // max neighbors of any site
  int *numneigh;               // # of neighbors of each owned site
  int **neighbor;              // list of neighbors of each owned site
                               // neighbor[i][j] =
                               //   local index of Jth neigh of Ith owned site
                               //   can index an owned or ghost site

  int nbasis;                  // basis atoms for regular lattices
  int ***cmap;                 // connectivity map for regular lattices

  struct Ghost {
    int id,proc,index;
    double x,y,z;
  };

  struct {                     // storage for a single site with general data
    int *ivalue;
    double *dvalue;
  } onesite;

  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition

  FILE *fp;
  double *dbuf;
  int maxdumpbuf;

  class RandomPark *random;
  class CommLattice *comm;

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}

  void options(int, char **);
  void create_lattice();

  void structured_lattice();
  void random_lattice();
  void file_lattice();
  void ghosts_from_connectivity();
  void ghosts_within_cutoff();

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
