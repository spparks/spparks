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

#ifndef APP_PORE2_H
#define APP_PORE2_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppPore2 : public AppLattice {
 public:
  AppPore2(class SPPARKS *, int, char **);
  ~AppPore2();
  void init_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int *sites;
  int *check;
  double *ecoord;
  double kboltz;

  void input_app(char *, int, char **);

  struct Event {           // one event for an owned site
    int partner;           // local ID of exchange partner
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, double);
};

}

#endif
