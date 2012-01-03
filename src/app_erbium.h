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
AppStyle(erbium,AppErbium)

#else

#ifndef SPK_APP_ERBIUM_H
#define SPK_APP_ERBIUM_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppErbium : public AppLattice {
  friend class DiagErbium;

 public:
  AppErbium(class SPPARKS *, int, char **);
  ~AppErbium();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle;
  int *type,*element;      // variables on each lattice site
  int firsttime;

  int *esites;
  int *echeck;

  int none,ntwo,nthree;
  double *srate,*drate,*trate;
  double *spropensity,*dpropensity,*tpropensity;
  int *stype,**dtype,**ttype;
  int *sinput,**dinput,**tinput;
  int *soutput,**doutput,**toutput;
  int *scount,*dcount,*tcount;

  struct Event {           // one event for an owned site
    int style;             // reaction style = SINGLE,DOUBLE,TRIPLE
    int which;             // which reaction of this type
    int jpartner,kpartner; // which J,K neighbors of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, int, double, int, int);
  void grow_reactions(int);
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

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Temperature cannot be 0.0 for app erbium

UNDOCUMENTED

*/
