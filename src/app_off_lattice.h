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

#ifndef SPK_APP_OFF_LATTICE_H
#define SPK_APP_OFF_LATTICE_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS_NS {

class AppOffLattice : public App {
  friend class CommOffLattice;

 public:
  int nglobal;                 // global # of sites
  int nlocal;                  // # of sites I own
  class CommOffLattice *comm;

  AppOffLattice(class SPPARKS *, int, char **);
  virtual ~AppOffLattice();
  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();
  void *extract(char *);

  // pure virtual functions, must be defined in child class

  virtual double site_energy(int) = 0;
  virtual void site_event_rejection(int, class RandomPark *) = 0;
  virtual double site_propensity(int) = 0;
  virtual void site_event(int, class RandomPark *) = 0;

  // virtual functions, may be overridden by child class

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}
  void virtual setup_app() {}

  enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
	 FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D,FILENAME};

 protected:
  int me,nprocs;
  int naccept,nattempt;       // number of accepted and attempted events
  int nsweeps;                // number of sweeps performed
  double temperature,t_inverse;  // temperature settings
  double dt_sweep;            // rKMC time for nglobal attemped events
  double dt_rkmc;             // rKMC time for one pass thru all sectors
  double dt_kmc;              // KMC time for one pass thru all sectors

  class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
  int *sitelist;               // randomized list of site indices

  int latstyle;               // lattice creation params
  double latconst;
  int dimension;
  int nx,ny,nz;
  int nrandom;
  double cutoff;
  char *latfile;
  char *infile;

  int px_user,py_user,pz_user;

  double delpropensity;        // distance away needed to compute propensity
  double delevent;             // distance away affected by an event
  int allow_kmc;               // 1 if app supports KMC
  int allow_rejection;         // 1 if app supports rejection KMC

  int sweepflag;               // set if rejection KMC solver
  int sectorflag;              // 1 if partition my domain into sectors
  int nsector;                 // 1,2,4,8 = # of sectors
  int nsector_user;            // 0 if default, else 2,4,8

  double Ladapt;               // adaptive sector time increments for KMC
  double tstop;                // requested time increment in sector
  double nstop;                // requested events per site in sector

  double xprd,yprd,zprd;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // simulation box bounds
  double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;    // my portion of box
  int nx_procs,ny_procs,nz_procs;   // procs in each dim of lattice partition

  int nghost;                  // # of ghost sites I store

                               // these arrays stored for owned + ghost sites
  int *id;                     // global ID (1-N) of site
  int *owner;                  // proc who owns the site
  double **xyz;                // coords of site

                               // per-site storage for owned + ghost sites
  int sitecustom;              // 0/1 for default or customized
  int ninteger,ndouble;        // # of int/double per site, 0,0 = just lattice
  int *site;                   // default = single int value
  int **iarray;                // one or more ints per site
  double **darray;             // one or more doubles per site

  double *propensity;          // probabilities for each owned site
  int *i2site;                 // mapping of owned lattice to site index

                               // neighbor list info
  int numneigh;
  int *neighs;

  int nbasis;                  // basis atoms for regular lattices

  struct Set {                 // subset of particles I own
    int nlocal;                // # of owned sites in set
    int nborder;               // # of sites with non-set site as neighbor
    int nselect;               // # of selections from set for rKMC
    int nloop;                 // # of loops over set for rKMC
    int *border;               // lattice index for each border site
    int *bsites;               // list of border sites to pass to solver
    class Solve *solve;        // KMC solver
    double *propensity;
    int *site2i;               // map from set sites to lattice index
    int *i2site;               // map from lattice index to set sites
  };
  Set *set;                    // list of subsets
  int nset;                    // # of subsets of lattice sites


  void iterate_kmc_global(double);
  void iterate_kmc_sector(double);
  void iterate_rejection(double);

  void neighbor(int, double) {}

  void options(int, char **);
  void create_domain();
  void structured_lattice();
  void random_lattice();
  void file_lattice();
  void read_file();

  void create_set(int, int);

  void set_sector(int, char **);
  void set_sweep(int, char **);
  void set_temperature(int, char **);

  void stats(char *);
  void stats_header(char *);
};

}

#endif
