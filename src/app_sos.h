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
AppStyle(sos,AppSOS)

#else

#ifndef SPK_APP_SOS_H
#define SPK_APP_SOS_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppSOS : public AppLattice {
 public:
  AppSOS(class SPPARKS *, int, char **);
  ~AppSOS();
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {};
  double site_propensity(int);
  void site_event(int, class RandomPark *);

  const int *get_lattice() { return iarray[0]; };
  int get_nlocal() { return nlocal; };
  double *inputlat;

 private:
  double boltz;
  int *height;
  int *sites;
  int *check;
  double stepheight,bondeng,fullp;
  double amp,xwl,zwl;
  int instyle;
  double *invalues;
  int firsttime;
  double tscale_inverse;

  struct Event {           // one event for an owned site
    double propensity;     // propensity of this event
    int partner;           // local ID of exchange partner
    int next;              // index of next event for this site
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, double, int);

  void create_height(double, double, double);
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

*/
