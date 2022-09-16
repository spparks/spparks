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

#ifdef APP_CLASS
AppStyle(diffusion,AppDiffusion)

#else

#ifndef SPK_APP_DIFFUSION_H
#define SPK_APP_DIFFUSION_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppDiffusion : public AppLattice {
  friend class DiagDiffusion;

 public:
  AppDiffusion(class SPPARKS *, int, char **);
  ~AppDiffusion();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

  void app_update(double);

 private:
  int engstyle,hopstyle;
  int allocated;
  int *esites,*psites;
  int *echeck,*pcheck;
  double *ecoord;

  int dimension;
  int *lattice;

  struct Event {           // one event for an owned site
    double propensity;     // propensity of this event
    int destination;       // local ID of destination site
    int style;             // nearest-neigh hop or Schwobel hop
    int next;              // index of next event for this site
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  int barrierflag;          // energy barriers on or off
  double **hbarrier;
  double **sbarrier;

  int nsmax,nsmin;          // Schwoebel hop params
  int *hopsite;             // list of possible hops for one site
  int *mark;                // flagged sites
  int *marklist;            // list of flagged sites

  // deposition parameters and data structs

  int depmode;             // deposition: DEP_NONE or DEP_EVENT or DEP_BATCH
  double deprate_total;    // deposition rate for entire system
  double deprate;          // deposition rate for one proc
  double thetalo,thetahi;  // deposition params
  double d0;
  int coordlo,coordhi;
  double dir[3];

  int nbatch;               // # of batch deposition events by this proc
  int maxbatch;             // max # of batch depositions data structs can hold
  class RandomPark *ranbatch;  // RNG for batch deposition (same on all procs)
  double **startpos;        // starting points for batch deposition
  int *elist;               // list of my sites eligible for deposition

  struct DepInfo {
    int proc;
    int site;
    double distance;
  };

  DepInfo *depinfo;        // allocated for batch depositions
  DepInfo *depinfo_copy;   // used to merge two sets of depinfo
  
  // stats

  int ndeposit,ndeposit_failed;
  int nfirst,nsecond;

  // methods

  double site_propensity_no_energy(int);
  double site_propensity_linear(int);
  double site_propensity_nonlinear(int);

  void site_event_linear(int, class RandomPark *);
  void site_event_nonlinear(int, class RandomPark *);
  void update_propensities(int, int);

  int neighbor2(int, int *);
  int neighbor3(int, int *);
  int neighbor4(int, int *);

  int ncoord(int);
  void clear_events(int);
  void add_event(int, int, double, int);

  int schwoebel_enumerate(int, int *);
  int find_deposition_site(class RandomPark *);
  int exceeds_limit(int, double *, double &);
  double distsq_to_line(int, double *, int, int, double &);
  void allocate_data();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot use %s command until sites exist

This command requires sites exist before using it in an input script.

E: Can only use ecoord command with app_style diffusion nonlinear

Self-explanatory.

E: Cannot define Schwoebel barrier without Schwoebel model

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Cannot perform deposition in parallel

UNDOCUMENTED

E: Cannot perform deposition with multiple sectors

UNDOCUMENTED

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Did not reach event propensity threshhold

UNDOCUMENTED

E: BAD DONE

UNDOCUMENTED

*/
