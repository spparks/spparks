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

#ifdef AppClass
AppStyle(diffusion,AppDiffusion)

#else

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppDiffusion : public AppLattice {
  friend class DiagDiffusion;

 public:
  AppDiffusion(class SPPARKS *, int, char **);
  ~AppDiffusion();
  void input_app(char *, int, char **);
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle,hopstyle,geomstyle;
  int *esites,*psites;
  int *echeck,*pcheck;
  double *ecoord;

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

  int depflag;             // deposition on or off
  double deprate,thetalo,thetahi;
  double d0;
  int coordlo,coordhi;
  double dir[3];

  int barrierflag;          // energy barriers on or off
  double **hbarrier;
  double **sbarrier;

  int nsmax,nsmin;          // Schwoebel hop params
  int *hopsite;             // list of possible hops for one site
  int *mark;                // flagged sites
  int *marklist;            // list of flagged sites

  int ndeposit,ndeposit_failed;  // stats
  int nfirst,nsecond;

  double site_propensity_no_energy(int);
  double site_propensity_linear(int);
  double site_propensity_nonlinear(int);

  void site_event_linear(int, class RandomPark *);
  void site_event_nonlinear(int, class RandomPark *);

  int ncoord(int);

  void clear_events(int);
  void add_event(int, int, double, int);

  int find_deposition_site(class RandomPark *);
  int exceed_limit(int, double *, double &);
  double distsq_to_line(int, double *, int, int, double &);
  void bounds(char *, int, int &, int &);
  int schwoebel_enumerate(int, int *);
};

}

#endif
