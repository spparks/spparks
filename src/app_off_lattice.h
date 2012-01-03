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
  class CommOffLattice *comm;

  AppOffLattice(class SPPARKS *, int, char **);
  virtual ~AppOffLattice();
  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();

  void grow(int);
  void add_site(tagint, double, double, double);
  void add_values(int, char **x);

  // pure virtual functions, must be defined in child class

  virtual void grow_app() = 0;
  virtual double site_energy(int) = 0;
  virtual void site_event_rejection(int, class RandomPark *) = 0;
  virtual double site_propensity(int) = 0;
  virtual void site_event(int, class RandomPark *) = 0;

  // virtual functions, may be overridden by child class

  void virtual input_app(char *, int, char **);
  void virtual init_app() {}
  void virtual setup_app() {}
  virtual void *extract_app(char *) {return NULL;}

  enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
	 FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D,FILENAME};

 protected:
  int me,nprocs;

  int dimension;              // domain settings
  double xprd,yprd,zprd;
  double subxlo,subylo,subzlo;
  double subxhi,subyhi,subzhi;

  bigint naccept,nattempt;    // number of accepted and attempted events
  int nsweeps;                // number of sweeps performed
  double temperature,t_inverse;  // temperature settings
  double dt_sweep;            // rKMC time for nglobal attemped events
  double dt_rkmc;             // rKMC time for one pass thru all sectors
  double dt_kmc;              // KMC time for one pass thru all sectors

  class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
  int *sitelist;               // randomized list of site indices

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

  int nmax;                    // max # of sites per-site arrays can store
  int size_one;                // quantities to comm per site

                               // arrays for owned + ghost sites
  int *bin;                    // bin the site is in
  int *next,*prev;             // ptrs to next,previous site in bin, -1 if none
  int *nextimage;              // index of next image of site, -1 if none

  int *site2i;                 // list of site indices that are in sector
  int *in_sector;              // 1 if site is in sector, 0 if not

  double *propensity;          // probabilities for each owned site

                               // neighbor list info
  int numneigh;
  int *neighs;
  int nstencil;                // # of surrounding bins
  int *stencil;                // offsets for surrounding bins

  struct Set {                 // subset of particles I own
    double xlo,xhi;            // sector bounds
    double ylo,yhi;
    double zlo,zhi;
  };
  Set *set;                    // list of subsets
  int nset;                    // # of subsets of sites
  int activeset;               // which set is currently being looped over

                               // bin info
  int nbins;                   // total # of bins on this proc
  int nbinx,nbiny,nbinz;       // # of bins in each dim on this proc w/ ghosts
  double binx,biny,binz;       // bin sizes
  double invbinx,invbiny,invbinz;  // inverse bin sizes

  int *binhead;                // index of 1st particle in bin, -1 if none
  int *binflag;                // 0 if interior bin, 1 if edge, 2 if ghost

                               // image bins are the ghost copies of
                               // an edge bin, can be several per edge bin
  int *nimages;                // # of image bins for each edge bin
  int **imageproc;             // proc that owns the edge's ghost image bin
  int **imageindex;            // index of ghost image bin on owning proc

                               // master bin is the unique owned bin
                               // that corresponds to a ghost bin
  int *ghostproc;              // proc that owns a ghost's master bin
  int *ghostindex;             // index of master bin on owning proc

  int **pbcoffset;             // periodic offsets of each ghost bin in 3 dims
                               // offset to add to original site -> ghost site
                               // is 0 if ghost bin is not across PBC

  int nfree;                   // # of free locations in ghost particle list
  int freehead;                // 1st free location, -1 if none

  void iterate_kmc_global(double);
  void iterate_kmc_sector(double);
  void iterate_rejection(double);
  void rkmc_params(int, int &, int &);

  void init_bins();
  void init_stencil();

  void neighbor(int, double);
  void move(int);

  int site2bin(int);
  void delete_from_bin(int, int);
  void add_to_bin(int, int);
  void delete_images(int);
  void add_images(int, int);
  void add_to_free(int);
  void delete_all_ghosts();
  int new_ghost_site();
  int new_owned_site();
  int delete_owned_site(int);
  void add_free(int);
  void add_image_bins(int, int, int, int);
  int neighproc(int, int, int, int);
  int inside_sector(int);

  void check(char *, int, int);
  void create_set(int, int);

  void set_sector(int, char **);
  void set_sweep(int, char **);
  void set_temperature(int, char **);

  void stats(char *);
  void stats_header(char *);
};

}

#endif

/* ERROR/WARNING messages:

E: App needs a KMC or rejection KMC solver

You must define either a solver or sweep option.

E: App cannot use both a KMC and rejection KMC solver

You cannot define both a solver and sweep option.

E: KMC events are not implemented in app

Not every application supports KMC solvers.

E: Rejection events are not implemented in app

Self-explanatory.

E: Cannot use KMC solver in parallel with no sectors

Self-explanatory.

E: Cannot use random rejection KMC in parallel with no sectors

Self-explanatory.

E: Cannot use raster rejection KMC in parallel with no sectors

Self-explanatory.

E: App did not set dt_sweep

Internal SPPARKS error.

E: Invalid number of sectors

Self-explanatory.

E: Choice of sector stop led to no rKMC events

Self-explanatory.

E: SITE MISMATCH

UNDOCUMENTED

E: GHOST IN OWNED BIN

UNDOCUMENTED

E: LINK MISMATCH

UNDOCUMENTED

E: BIN MISMATCH

UNDOCUMENTED

E: SITES NOT IN BINS

UNDOCUMENTED

E: COUNT MISMATCH

UNDOCUMENTED

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Application cutoff is too big for processor sub-domain

There must be at least 2 bins per processor in each dimension
where sectoring occurs.

E: BAD STENCIL

UNDOCUMENTED

E: Too many neighbors per site

Internal SPPARKS error.

E: PBC remap of site failed

Internal SPPARKS error.

E: Site not in my bin domain

Internal SPPARKS error.

E: Adding site to illegal bin

Internal SPPARKS error.

E: Adding site to bin it is not in

Internal SPPARKS error.

E: Per-processor system is too big

UNDOCUMENTED

*/
