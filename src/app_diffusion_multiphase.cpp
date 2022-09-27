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
#include "string.h"
#include "stdlib.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "app_diffusion_multiphase.h"

using namespace SPPARKS_NS;
using std::map;
using std::set;

enum{LINEAR};

#define DELTAEVENT 100000

// This app is based on the diffusion app (for Kawasaki dynamics)
// These are the significant changes and simplifications
// 1. No Schowebel hops, only exchanges to site neighbors
// 2. Allow for an arbitrary number of phases (species for an atomic model)
// 3. Only the linear energy style is supported
// 4. One or more phases can be pinned, meaning they do not diffuse

/* ---------------------------------------------------------------------- */

AppDiffusionMultiphase::AppDiffusionMultiphase(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg), phase_labels(), is_pinned(), weights()
{
  // need to double check these values

  ninteger = 1;
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 0;
  numrandom = 1;

  // no args for this app

  if (narg > 1) error->all(FLERR,"Illegal app_style command");

  engstyle = LINEAR;

  create_arrays();
  esites = NULL;
  echeck = NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;

  allocated = 0;
}

/* ---------------------------------------------------------------------- */

AppDiffusionMultiphase::~AppDiffusionMultiphase()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

  if (strcmp(command,"diffusion/multiphase") == 0) 
    parse_diffmultiphase(narg,arg);
  else error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphase::parse_diffmultiphase(int narg, char **arg)
{
   // 2 args: diffusion/multiphase phase <int value>
   // 2 args: diffusion/multiphase pin <int value>
   // 5 args: diffusion/multiphase weight <double> pair <int,int>

   if (narg < 2) 
     error->all(FLERR,"Illegal diffusion/multiphase command");

   if (strcmp(arg[0],"phase")==0){
      if (narg != 2) 
        error->all(FLERR,"Illegal diffusion/multiphase phase command: "
                   "num args != 2");
      int phase = std::atoi(arg[1]);
      phase_labels.insert(phase);
      // phases are not pinned by default
      is_pinned[phase] = false;
      if (phase < 1) 
        error->all(FLERR,"Illegal diffusion/multiphase phase value: "
                   "must be >= 1");
   } else if (strcmp(arg[0],"pin") == 0) {
      if (narg != 2) 
        error->all(FLERR,"Illegal diffusion/multiphase pin command: "
                   "num args != 2");
      int pin_phase = std::atoi(arg[1]);
      phase_labels.insert(pin_phase);
      is_pinned[pin_phase] = true;
      if (pin_phase < 1)
        error->all(FLERR,"Illegal diffusion/multiphase pin value: must be >= 1");
   } else if (strcmp(arg[0],"weight")==0){
      if (narg != 5) 
        error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                   "num args != 5");
      double w = std::atof(arg[1]);
      if (strcmp(arg[2],"pair") == 0) {
         int p1 = std::atoi(arg[3]);
         int p2 = std::atoi(arg[4]);
         if (p1 == p2) 
           error->all(FLERR,"Cannot set diffusion/multiphase weight for "
                      "pair of identical phases");
         if (p1 < 1 || p2 < 1) 
           error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                      "phases must be >= 1");
         if (w < 0.0) 
           error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                      "weight must be >= 0");
         weights[{p1,p2}] = w;
         weights[{p2,p1}] = w;
      } else 
        error->all(FLERR,"Illegal diffusion/multiphase weight command: "
                   "expected keyword pair");
   } else error->all(FLERR,"Illegal diffusion/multiphase command: "
                     "expected phase, pin, or weight");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::grow_app()
{
  lattice = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::init_app()
{
   if (!allocated) allocate_data();
   allocated = 1;

   dimension = domain->dimension;
   dt_sweep = 1.0/maxneigh;

   {
     // insure all site values are in set of phase labels, otherwise error

     std::set<int>::iterator not_found=phase_labels.end();
     int flag = 0;
     for (int i = 0; i < nlocal; i++)
       if (not_found == phase_labels.find(lattice[i])) flag=1;
     int flagall;
     MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
     if (flagall) error->all(FLERR,"One or more sites have invalid values");
   }

   {
     // create default weights if they do not already exist
     // only unlike phases contribute to energy

     std::map<std::pair<int,int>,double>::iterator not_found=weights.end();
     for (auto p : phase_labels) {
       for (auto q : phase_labels) {
         if (p==q) continue;
         std::map<std::pair<int,int>,double>::iterator i = weights.find({p,q});
         // if pair not found, set default weight = 1.0
         if (not_found==i) {
           weights[{p,q}] = 1.0;
           weights[{q,p}] = 1.0;
         }
       }
     }
   }
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusionMultiphase::site_energy(int i)
{
  // energy of site = linear sum of bond weights
  // only unlike phase pairs contribute

   double energy(0.0);
   int ip = lattice[i];
   for (int j = 0; j < numneigh[i]; j++){
      int nj = neighbor[i][j];
      int jp = lattice[nj];
      if (ip == jp) continue;
      energy += weights[{ip,jp}];
   }

   // each site carries half the interaction energy
   // neighbor sites carry the other half 
   return 0.5*energy;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   null bin extends to size maxneigh
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event_rejection(int i, RandomPark *random)
{
  double einitial,edelta;
  int i_old, j_old;

  // pinned sites can't exchange

  if (is_pinned[lattice[i]]) return;

  // need to double check neighborhood

  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];

  // if site j pinned or if site i and site j have same spin, just return

  if (is_pinned[lattice[j]] || lattice[j] == lattice[i]) return;

  i_old = lattice[i];
  j_old = lattice[j];
  
  // accept or reject via energy model

  int hop = 0;
  einitial = site_energy(i)+site_energy(j);

  lattice[i] = j_old;
  lattice[j] = i_old;

  // compute energy difference from exchange

  edelta = site_energy(i)+site_energy(j) - einitial;

  // if edelta is negative, accept
  // otherwise if temperature is non-zero, can still accept

  if (edelta <= 0.0) hop = 1;
  else if (temperature > 0.0) {
    if (random->uniform() < exp(-1.0*edelta*t_inverse)) hop = 1;
  }
    
  if (hop) {
    naccept++;
  } else {
    lattice[i] = i_old;
    lattice[j] = j_old;
  }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusionMultiphase::site_propensity(int i)
{
  return site_propensity_linear(i);
}

/* ---------------------------------------------------------------------- */

double AppDiffusionMultiphase::site_propensity_linear(int i)
{
  int j,k, i_old, j_old;
  double e0,einitial,edelta,probone,proball;

  // add event if neighbors are dissimilar and not pinned

  clear_events(i);

  if (is_pinned[lattice[i]]) return 0.0;
  
  i_old = lattice[i];

  // loop over all possible hops, go through neighbor shell

  // einitial = site_energy(i);
  e0 = site_energy(i);
  proball = 0.0;
  probone = 0.0;
  
  // this is similar to the Potts approach

  for (k = 0; k < numneigh[i]; k++) {
    j = neighbor[i][k];
    j_old = lattice[j];
    
    if (lattice[j] == lattice[i] || is_pinned[lattice[j]]) continue;

    //einitial = site_energy(i_old)+site_energy(j_old);
    einitial = e0+site_energy(j);
    
    // exchange values and check energy

    lattice[i] = j_old;
    lattice[j] = i_old;
    edelta = site_energy(i)+site_energy(j) - einitial;
    
    // if energy is non-zero, add as possible event

    if (edelta <= 0.0) probone = 1.0;
    else if (temperature > 0.0) {
      probone = exp(-1.0*edelta*t_inverse);
    }
    
    if (probone > 0.0) {
      add_event(i,j,probone);
      proball += probone;
    }

    lattice[i] = i_old;
    lattice[j] = j_old;
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event(int i, class RandomPark *random)
{
  return site_event_linear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusionMultiphase::site_event_linear(int i, class RandomPark *random)
{
  int j,k,m,isite,i_old,j_old;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one(FLERR,"Did not reach event propensity threshhold");
  }

  // exchange values

  j = events[ievent].destination;
  i_old = lattice[i];
  j_old = lattice[j];
  lattice[i] = j_old;
  lattice[j] = i_old;

  // compute propensity changes for self and swap site and their neighs
  // 1,2 neighs for NNHOP and 1,2,3 neighs for SCHWOEBEL
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
  
  // update site i's neighbors, this will include the exchanged site

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    // not quite sure what this does
    if (isite < 0) continue;
    // add to update list
    esites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
    echeck[isite] = 1;
  }
  
  // update exchanged site's neighbors
  // avoid any that have already been found

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    // not quite sure what this does
    if (isite < 0) continue;
    // make sure site is not already updated
    if (echeck[isite] == 1) continue;
    // add to update list
    esites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
    echeck[isite] = 1;
  }
  
  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (k = 0; k < nsites; k++) echeck[esites[k]] = 0;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::clear_events(int i)
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

void AppDiffusionMultiphase::add_event(int i, int destination, 
			      double propensity)
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
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppDiffusionMultiphase::allocate_data()
{
  // for linear:
  //   make esites large enough for 1 sites and their 1,2 neighbors

  if (engstyle == LINEAR) {
    int emax = 1 + maxneigh*2;
    esites = new int[2*emax];
  }

  echeck = new int[nlocal+nghost];

  memory->create(firstevent,nlocal,"app:firstevent");
}
