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

/* ----------------------------------------------------------------------
   Contributing authors: Greg Wagner (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_sos.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED};
enum{XSIN,STRING,NONE};

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppSOS::AppSOS(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  create_arrays();

  // parse arguments

  if (narg < 2) error->all(FLERR,"Illegal app_style command");
  
  bondeng = atof(arg[1]);

  instyle = NONE;
  invalues = NULL;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"xsin") == 0) {
      if (iarg+4 < narg) error->all(FLERR,"Illegal app_style command");
      instyle = XSIN;
      amp = atof(arg[iarg+1]);
      xwl = atof(arg[iarg+2]);
      zwl = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"string") == 0) {
      if (iarg+2 < narg) error->all(FLERR,"Illegal app_style command");
      instyle = STRING;
      invalues = (double *) strtol(arg[iarg+1],NULL,16);
      iarg += 2;
    } else error->all(FLERR,"Illegal app_style command");
  }

  // defaults

  boltz = 1.0;
  stepheight = 1.0;

  firsttime = 1;
  sites = check = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;
}

/* ---------------------------------------------------------------------- */

AppSOS::~AppSOS()
{
  memory->sfree(invalues);

  delete [] sites;
  delete [] check;
  memory->sfree(events);
  memory->sfree(firstevent);
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppSOS::grow_app()
{
  height = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   no check for site validity, since any int value is OK
------------------------------------------------------------------------- */

void AppSOS::init_app()
{
  // NOTE: probably remove STRING option
  // NOTE: how to insure invalues is right length?

  if (firsttime) {
    firsttime = 0;

    if (instyle == XSIN) create_height(xwl,zwl,amp);
    else if (instyle == STRING)
      for (int i = 0; i < nlocal; i++)
	height[i] = (int) invalues[i];

    check = new int[nlocal+nghost];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");

    // sites must be large enough for 2 sites and their 1,2 nearest neighbors

    int nmax = 1 + maxneigh + maxneigh*maxneigh;
    sites = new int[2*nmax];
  }

  // prefactor on propensity is 1/N_neighbors

  fullp = 1.0/maxneigh;

  tscale_inverse = t_inverse/boltz;
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppSOS::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) check[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppSOS::site_energy(int i)
{
  int j,isite,jsite;

  isite = height[i];

  int eng = 0;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    jsite = height[j];
    if (j == i) continue;
    eng = eng + abs(isite-jsite);
  }

  return (double) 0.5*eng*bondeng*stepheight; 
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppSOS::site_propensity(int i)
{
  int mystate,neighstate,ineigh,j;
  double einitial,efinal,einitiali,proball,probone;

  mystate = height[i];

  clear_events(i);

  proball = 0.0;
  einitiali = 2.0*site_energy(i);
  for (ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    neighstate = height[j];
    if (j == i) continue;
    einitial = einitiali + 2.0*site_energy(j)-
               bondeng*abs(mystate-neighstate)*stepheight;
    height[i] = mystate-1;
    height[j] = neighstate+1;
    efinal = 2.0*(site_energy(i)+site_energy(j))-
             bondeng*abs(height[i]-height[j])*stepheight;

    if (efinal <= einitial) probone = fullp;
    else probone = fullp*exp((einitial-efinal)*tscale_inverse);
    if (probone > 0.0) {
      add_event(i,j,probone,0);
      proball += probone;
    }
    height[i] = mystate;
    height[j] = neighstate;
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppSOS::site_event(int i, class RandomPark *random)
{
  int j,k,kk,m,mm,isite,nsites;

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
  }
  j = events[ievent].partner;
  height[i] = height[i]-1;
  height[j] = height[j]+1;

  // compute propensity changes for self and swap site and their 1,2 neighs
  // ignore update of sites with isite < 0
  // use check[] to avoid resetting propensity of same site

  nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  sites[nsites++] = isite;
  check[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
    sites[nsites++] = isite;
    check[isite] = 1;
  }

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && check[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      check[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && check[isite] == 0) {
        propensity[isite] = site_propensity(mm);
        sites[nsites++] = isite;
        check[isite] = 1;
      }
    }
  }

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    if (isite >= 0 && check[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      check[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && check[isite] == 0) {
        propensity[isite] = site_propensity(mm);
        sites[nsites++] = isite;
        check[isite] = 1;
      }
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

void AppSOS::clear_events(int i)
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
   event = exchange with partner with probability = propensity
------------------------------------------------------------------------- */

void AppSOS::add_event(int i, int partner, double propensity, int incre)
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
  events[freeevent].partner = partner;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
  create height values
------------------------------------------------------------------------- */
 
void AppSOS::create_height(double xwl, double zwl, double amp)
{
  double x,y,z;
  double ylocalx,ylocalz,ylocal;

  for (int i = 0; i < nlocal; i++) {
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    ylocalx = amp*sin(2.0*M_PI*x/xwl);
    if (zwl < 1.0e10) ylocalz = amp*sin(2.0*M_PI*z/zwl);
    else ylocalz = amp*sin(2.0*M_PI*x/xwl);
    if (ylocalx < ylocalz) ylocal = ylocalx;
    else ylocal = ylocalz;
    if (ylocal > 0.0) ylocal += 0.5;
    else if (ylocal < 0.0) ylocal -= 0.5;
    height[i] = (int) ylocal;
  }
}
