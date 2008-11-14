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
  // esites must be large enough for 2 sites and 1st/2nd nearest neighbors

  create_lattice();
  esites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  echeck = NULL;

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
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
}

/* ---------------------------------------------------------------------- */

void AppErbium::init_app()
{
  delete [] echeck;
  echeck = new int[nlocal];
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

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
   propensity for one event is based on einitial,efinal
------------------------------------------------------------------------- */

double AppErbium::site_propensity(int i)
{
  int j;

  // possible events from tabulated list

  clear_events(i);

  if (element[i] != ERBIUM) return 0.0;   // not correct, should allow all sites to have events, ERBIUM only handles 2-hop exchange events

  double einitial,efinal,probone;
  double proball = 0.0;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANCY) {
      einitial = site_energy(i) + site_energy(j);

//      lattice[i] = VACANT;
//      lattice[j] = OCCUPIED;
      efinal = site_energy(i) + site_energy(j);

      if (efinal <= einitial) probone = 1.0;
      else if (temperature > 0.0) probone = exp((einitial-efinal)*t_inverse);
      else probone = 0.0;

      if (probone > 0.0) {

	// also put event ID into event list?

	add_event(i,j,probone);
	proball += probone;
      }

      //      lattice[i] = OCCUPIED;
      // lattice[j] = VACANT;
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
  int j,jj,k,kk,m,mm;

  // pick one event from total propensity for this site
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  j = events[ievent].partner;
  //  lattice[i] = VACANT;
  //lattice[j] = OCCUPIED;

  // compute propensity changes for self and swap site and their 1,2 neighs
  // use echeck[] to avoid resetting propensity of same site
  // loop over neighbors of site even if site itself is out of sector

  int nsites = 0;
  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
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

void AppErbium::add_event(int i, int partner, double propensity)
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
  events[freeevent].partner = partner;
  events[freeevent].next = firstevent[i];
  events[freeevent].propensity = propensity;
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}
