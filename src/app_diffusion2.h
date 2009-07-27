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

#ifndef APP_DIFFUSION2_H
#define APP_DIFFUSION2_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppDiffusion2 : public AppLattice {
  friend class DiagDeposition;

 public:
  AppDiffusion2(class SPPARKS *, int, char **);
  ~AppDiffusion2();
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
    int destination;       // local ID of destination site
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  int depflag;             // deposition on or off
  int ndeposit;
  int ndeposit_failed;
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

  double site_propensity_linear(int);
  double site_propensity_table(int);
  double site_propensity_nonlinear(int);

  void site_event_linear(int, class RandomPark *);
  void site_event_table(int, class RandomPark *);
  void site_event_nonlinear(int, class RandomPark *);

  int ncoord(int);

  void clear_events(int);
  void add_event(int, int, double);

  int find_deposition_site(class RandomPark *);
  int exceed_limit(int, double *, double &);
  double distsq_to_dir(int, double *, int, int, double &);
  void bounds(char *, int, int &, int &);
  int schwoebel_enumerate(int, int *);
};

}

#endif
