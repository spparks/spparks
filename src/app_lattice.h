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

#ifndef SPK_APP_LATTICE_H
#define SPK_APP_LATTICE_H

#include "stdio.h"
#include "app.h"

#include <stack>

namespace SPPARKS_NS {

class AppLattice : public App {
  friend class Dump;
  friend class CommLattice;
  friend class DiagEnergy;
  friend class DiagErbium;
  friend class DiagCluster;

 public:
  int nglobal;                 // global # of sites
  int nlocal;                  // # of sites I own
  class CommLattice *comm;

  AppLattice(class SPPARKS *, int, char **);
  virtual ~AppLattice();
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
  virtual void *extract_app(char *) {return NULL;}

  virtual void push_new_site(int, int *, int, std::stack<int>*);
  virtual void push_connected_neighbors(int, int *, int, std::stack<int>*);
  virtual void connected_ghosts(int, int *, class Cluster *, int);

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

  int latstyle;               // lattice creation params
  double latconst;
  int dimension;
  int nx,ny,nz;
  int nrandom;
  double cutoff;
  char *latfile;
  char *infile;
  int px_user,py_user,pz_user;

  int delpropensity;           // # of sites away needed to compute propensity
  int delevent;                // # of sites away affected by an event
  int numrandom;               // # of RN used by rejection routine
  int allow_kmc;               // 1 if app supports KMC
  int allow_rejection;         // 1 if app supports rejection KMC
  int allow_masking;           // 1 if app supports rKMC masking

  int sweepflag;               // 1 if rejection KMC solver
  int sectorflag;              // 1 if partition my domain into sectors
  int nsector;                 // 1,2,4,8 = # of sectors
  int nsector_user;            // 0 if default, else 2,4,8
  int ncolors;                 // # of colors, depends on lattice
  int bothflag;                // 1 if both sectors and colors

  class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
  class RandomPark *ranstrict; // RN generator for per-site strict rKMC
  int *siteseeds;              // per-site seeds for ransite
  int *sitelist;               // randomized list of site indices

  bool Lmask;                  // masking on/off
  char *mask;                  // size of nlocal + nghost sites

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
                               // cmap[nbasis][maxneigh][4]
                               // 0,1,2 = x,y,z offsets in unit cell
                               // 3 = which atom in offset unit cell

  struct Set {                 // subset of lattice sites I own
    int nlocal;                // # of owned sites in set
    int nselect;               // # of selections from set for rKMC
    int nloop;                 // # of loops over set for rKMC
    int nborder;               // # of sites with non-set site as neighbor
    int *border;               // lattice index for each border site
    int *bsites;               // list of border sites to pass to solver
    class Solve *solve;        // KMC solver
    double *propensity;        // propensities for set sites
    int *site2i;               // map from set sites to lattice index
    int *i2site;               // map from lattice index to set sites
  };
  Set *set;                    // list of subsets
  int nset;                    // # of subsets of lattice sites

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  void iterate_kmc_global(double);
  void iterate_kmc_sector(double);
  void iterate_rejection(double);

  typedef void (AppLattice::*FnPtrSweep)(int, int *);
  FnPtrSweep sweep;                         // ptr to sweep functions
  void sweep_nomask_nostrict(int, int *);
  void sweep_mask_nostrict(int, int *);
  void sweep_nomask_strict(int, int *);
  void sweep_mask_strict(int, int *);

  void options(int, char **);
  void create_lattice();
  void structured_lattice();
  void random_lattice();
  void file_lattice();
  void read_file();

  void ghosts_from_connectivity();
  void connectivity_within_cutoff();

  void create_set(int, int, int);
  int id2color(int);
  int find_border_sites(int);
  void boundary_clear_mask(int);

  void stats(char *);
  void stats_header(char *);

  void set_sector(int, char **);
  void set_sweep(int, char **);
  void set_temperature(int, char **);

  int connect(int, int);
  void offsets(double **);
  void offsets_2d(int, double **, double, double, double, double, int, int **);
  void offsets_3d(int, double **, double, double, double, double, double,
		  int, int **);
};

}

#endif
