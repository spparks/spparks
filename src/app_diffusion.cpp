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
#include "app_diffusion.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppDiffusion::AppDiffusion(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;

  // parse arguments

  if (narg < 3) error->all("Illegal app_style command");

  double fraction = atof(arg[1]);
  int seed = atoi(arg[2]);
  random = new RandomPark(seed);

  options(narg-3,&arg[3]);

  // define lattice and partition it across processors
  // sites must be large enough for 2 sites and their 1st/2nd nearest neighbors

  create_lattice();
  sites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  check = NULL;

  // event list

  events = NULL;
  firstevent = NULL;

  // initialize my portion of lattice
  // each site = 1 (vacancy) or 2 (occupied) with fraction occupied
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  std::map<int,int> hash;
  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));
  std::map<int,int>::iterator loc;

  int isite;
  for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
    if (random->uniform() < fraction) isite = 2;
    else isite = 1;
    loc = hash.find(iglobal);
    if (loc != hash.end()) lattice[loc->second] = isite;
  }
}

/* ---------------------------------------------------------------------- */

AppDiffusion::~AppDiffusion()
{
  delete random;
  delete [] sites;
  memory->sfree(check);
  memory->sfree(events);
  memory->sfree(firstevent);
}

/* ---------------------------------------------------------------------- */

void AppDiffusion::init_app()
{
  memory->sfree(check);
  check = (int *) memory->smalloc(nlocal*sizeof(int),"app:check");
  for (int i = 0; i < nlocal; i++) check[i] = 0;

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

double AppDiffusion::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != lattice[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
------------------------------------------------------------------------- */

void AppDiffusion::site_event_rejection(int i, RandomPark *random)
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

double AppDiffusion::site_propensity(int i)
{
  int j;

  // possible events = exchange with neighboring site different than self

  clear_events(i);

  int mystate = lattice[i];
  int neighstate;
  double proball = 0.0;
  double einitial,efinal,probone;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    neighstate = lattice[j];
    if (neighstate != mystate) {
      einitial = site_energy(i) + site_energy(j);
      lattice[i] = neighstate;
      lattice[j] = mystate;
      efinal = site_energy(i) + site_energy(j);
      if (efinal <= einitial) probone = 1.0;
      else if (temperature > 0.0) probone = exp((einitial-efinal)*t_inverse);
      if (probone > 0.0) {
	add_event(i,j,probone);
	proball += probone;
      }
      lattice[i] = mystate;
      lattice[j] = neighstate;
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   ignore neighbor sites that should not be updated (isite < 0)
------------------------------------------------------------------------- */

void AppDiffusion::site_event(int i, class RandomPark *random)
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
  int tmp = lattice[i];
  lattice[i] = lattice[j];
  lattice[j] = tmp;

  // compute propensity changes for self and swap site
  // 2nd neighbors of I,J could change their propensity
  // use check[] to avoid resetting propensity of same site
  // NOTE: should I loop over 2nd neighbors even if site itself is skipped?
  //   not if skipped b/c already seen, but maybe if out of sector

  int nsites = 0;
  int isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);
  check[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(j);
    check[isite] = 1;
  }

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite < 0) continue;
    if (check[isite] == 0) {
      sites[nsites++] = isite;
      propensity[isite] = site_propensity(m);
      check[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite < 0 || check[isite]) continue;
      sites[nsites++] = isite;
      propensity[isite] = site_propensity(mm);
      check[isite] = 1;
    }
  }

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    if (isite < 0) continue;
    if (check[isite] == 0) {
      sites[nsites++] = isite;
      propensity[isite] = site_propensity(m);
      check[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite < 0 || check[isite]) continue;
      sites[nsites++] = isite;
      propensity[isite] = site_propensity(mm);
      check[isite] = 1;
    }
  }

  solve->update(nsites,sites,propensity);

  // clear check array

  for (m = 0; m < nsites; m++) check[sites[m]] = 0;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusion::clear_events(int i)
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

void AppDiffusion::add_event(int i, int partner, double propensity)
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
