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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_erbium.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{NOOP,FCC,TETRA,OCTA};
enum{ZERO,ERBIUM,HYDROGEN,HELIUM,VACANCY};
#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppErbium::AppErbium(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 1;

  // parse arguments

  if (narg < 3) error->all("Illegal app_style command");

  double fraction = atof(arg[1]);
  int seed = atoi(arg[2]);
  random = new RandomPark(seed);

  options(narg-3,&arg[3]);

  // define lattice and partition it across processors
  // sites must be large enough for self and 1st neighbors

  create_lattice();
  sites = new int[1 + maxneigh];

  // assign variable names

  if (ninteger != 2 || ndouble != 0)
    error->all("Invalid site specification in app_style erbium");

  type = iarray[0];
  element = iarray[1];

  // event list

  events = NULL;
  firstevent = NULL;

  // initialize my portion of lattice
  // FCC sites are ERBIUM
  // OCTA sites are VACANCY
  // random fraction of TETRA sites are HYDROGEN
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  if (infile) read_file();

  else {
    std::map<int,int> hash;
    for (int i = 0; i < nlocal; i++)
      hash.insert(std::pair<int,int> (id[i],i));
    std::map<int,int>::iterator loc;
    
    int itype,flag;
    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      if (random->uniform() < fraction) flag = HYDROGEN;
      else flag = VACANCY;
      loc = hash.find(iglobal);
      if (loc != hash.end()) {
	itype = type[loc->second];
	if (itype == FCC) element[loc->second] = ERBIUM;
	else if (itype == OCTA) element[loc->second] = VACANCY;
	else if (itype == TETRA) element[loc->second] = flag;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

AppErbium::~AppErbium()
{
  delete random;
  delete [] sites;
  memory->sfree(events);
  memory->sfree(firstevent);
}

/* ---------------------------------------------------------------------- */

void AppErbium::init_app()
{
  memory->sfree(events);
  memory->sfree(firstevent);

  events = NULL;
  nevents = maxevent = 0;
  firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppErbium::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
------------------------------------------------------------------------- */

void AppErbium::site_event_rejection(int i, RandomPark *random)
{
  // event = exchange with random neighbor

  int iran = (int) (numneigh[i]*random->uniform());
  if (iran >= numneigh[i]) iran = numneigh[i] - 1;
  int j = neighbor[i][iran];

  double einitial = site_energy(i) + site_energy(j);

  int mystate = lattice[i];
  int neighstate = lattice[j];
  lattice[i] = neighstate;
  lattice[j] = mystate;

  double efinal = site_energy(i) + site_energy(j);

  // accept or reject the event

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = mystate;
    lattice[j] = neighstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = mystate;
    lattice[j] = neighstate;
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
   propensity for each possible event is input by user
------------------------------------------------------------------------- */

double AppErbium::site_propensity(int i)
{
  int j,k,m;

  // possible events from tabulated list

  clear_events(i);

  double proball = 0.0;

  // single-site events

  for (m = 0; m < nsingle; m++) {
    if (type[i] != stype[m] || element[i] != sinput[m]) continue;
    add_event(i,1,m,spropensity[m],-1,-1);
    proball += propensity[m];
  }

  // double-site events

  for (int j = 0; j < numneigh[i]; j++)
    for (m = 0; m < ndouble; m++) {
      if (type[i] != dtype[m][0] || element[i] != dinput[m][0]) continue;
      if (type[j] != dtype[m][1] || element[j] != dinput[m][1]) continue;
      add_event(i,2,m,dpropensity[m],j,-1);
      proball += dpropensity[m];
    }

  // triple-site events

  for (int j = 0; j < numneigh[i]; j++)
    for (int k = 0; k < numneigh[i]; k++) {
      if (j == k) continue;
      for (m = 0; m < ntriple; m++) {
	if (type[i] != ttype[m][0] || element[i] != tinput[m][0]) continue;
	if (type[j] != ttype[m][1] || element[j] != tinput[m][1]) continue;
	if (type[k] != ttype[m][2] || element[k] != tinput[m][2]) continue;
	add_event(i,3,m,tpropensity[m],j,k);
	proball += tpropensity[m];
      }
    }

  return proball;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   ignore neighbor sites that should not be updated (isite < 0)
------------------------------------------------------------------------- */

void AppErbium::site_event(int i, class RandomPark *random)
{
  int j,k,m;

  // pick one event from total propensity for this site
  // compare prob to threshhold, break when reach it to select event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform single, double, or triple event

  int type = events[ievent].type;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;

  if (type == 1) {
    element[i] = soutput[which];
  } else if (type == 2) {
    element[i] = doutput[which][0];
    element[neighbor[i][j]] = doutput[which][1];
  } else {
    element[i] = doutput[which][0];
    element[neighbor[i][j]] = toutput[which][1];
    element[neighbor[i][k]] = toutput[which][2];
  }

  // compute propensity changes for self and my first neighbors

  int nsites = 0;
  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  sites[nsites++] = isite;

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite >= 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
    }
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppErbium::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppErbium::add_event(int i, int type, int which, double propensity,
			  int jpartner, int kpartner)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].type = type;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}
