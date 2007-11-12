/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_LATTICE_H
#define APP_LATTICE_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppLattice : public App {
  friend class SweepLattice;
  friend class CommLattice;

 public:
  AppLattice(class SPK *, int, char **);
  virtual ~AppLattice();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  double virtual site_energy(int) = 0;
  int virtual site_pick_random(int,double) = 0;
  int virtual site_pick_local(int, double) = 0;
  double virtual site_propensity(int, int) = 0;
  void virtual site_event(int, int) = 0;
  void virtual site_clear_mask(char *, int) = 0;

 protected:
  int me,nprocs;
  int ntimestep,seed;
  int dump_style;
  double time,stoptime;
  double stats_time,stats_delta;
  double dump_time,dump_delta;
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

  int *lattice;                // lattice values for owned + ghost sites
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
  void stats();
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
