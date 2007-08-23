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

  int nglobal;
  int nlocal;
  int nghost;
  
  double xprd,yprd,zprd;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // lattice (box) bounds
  double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;    // my portion of box

  int *lattice;                // owned lattice + ghost lattice
  double *propensity;          // probability for each owned site

  int *id;
  int *owner;
  int *index;
  double **xyz;

  int *numneigh;
  int **neighbor;
  int maxconnect;
  int nbasis;
  int ***cmap;

  struct GhostRequest {
    int id,proc,index;
    double x,y,z;
  };

  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition
  int procwest,proceast;            // my neighbor procs
  int procsouth,procnorth;
  int procdown,procup;

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
  void structured_ghost();

  void random_lattice();
  void random_ghost();

  void file_lattice();

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
