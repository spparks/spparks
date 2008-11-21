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
#include "app_surface.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ZERO,FIXED,VACANT,OCCUPIED};
#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppSurface::AppSurface(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 2;
  delpropensity = 3;
  allow_metropolis = 0;

  boltz = 8.61734315e-5;
  vibrafreq=1.0e12;
  ebarrier=1.0;
  eSchwoebel=1.0e10;
  nminRegular = 1;
  nmaxSchwoebel = 3;
  nminSchwoebel = 2;
  bondener=-1.0;

  double xbottom,xheight,xmin,xmax,zbottom,zheight,zmin,zmax,yfix;
  int seed;

  xbottom=10000.0;
  xheight=0.0;
  xmin=0.0;
  xmax=100.0;
  zbottom=10000.0;
  zheight=0.0;
  zmin=0.0;
  zmax=100.0;
  yfix=-10000.0;
  seed=218392;

  // parse arguments

  int success = 0;

  for (int i = 1; i < narg; i++) {
     if (strcmp(arg[i],"xcos") == 0) {
        if (narg < i+5)  error->all("Illegal app_style command");
        else {
           xbottom = atof(arg[i+1]);
           xheight = atof(arg[i+2]);
           xmin = atof(arg[i+3]);
           xmax = atof(arg[i+4]);
           i = i+4;
        }
     }
     if (strcmp(arg[i],"zcos") == 0) {
        if (narg < i+5)  error->all("Illegal app_style command");
        else {
           zbottom = atof(arg[i+1]);
           zheight = atof(arg[i+2]);
           zmin = atof(arg[i+3]);
           zmax = atof(arg[i+4]);
           i = i+4;
        }
     }
     if (strcmp(arg[i],"yfix") == 0) {
        if (narg < i+2)  error->all("Illegal app_style command");
        else {
           yfix = atof(arg[i+1]);
           i = i+1;
        }
     }
     if (strcmp(arg[i],"seed") == 0) {
        if (narg < i+2)  error->all("Illegal app_style command");
        else {
           seed = atoi(arg[i+1]);
           i = i+1;
        }
     }
     if (strcmp(arg[i],"lattice") == 0) {
        if (narg < i+5)  error->all("Illegal app_style command");
        else {
           success = i;           
           i = i+4;
        }
     }
  }

  if (success == 0) error->all("Illegal app_style command");

  random = new RandomPark(seed);

  options(narg-success,&arg[success]);

  // define lattice and partition it across processors
  // sites must be large enough for 2 sites and their 1,2,3 nearest neighbors

  create_lattice();
  int nmax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
  sites = new int[2*nmax];
  sitesSchwoebel = new int[2*nmax];
  check = NULL;
  checkSchwoebel = NULL;

  // event list

  events = NULL;
  firstevent = NULL;

  // energy as a function of coordination

  ecoord = new double[maxneigh+1];
  for (int i = 0; i <= maxneigh; i++) ecoord[i] = 0.0;

  // initialize my portion of lattice
  // each site = 1 (vacancy) or 2 (occupied)
  // surface geometry defines occupied vs unoccupied

  double x,y,z;
  int isite;
  double ylocalx,ylocalz,ylocal;
  for (int i = 0; i < nlocal; i++) {
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    if (x <= xmin || x >= xmax) ylocalx = xbottom;
    else ylocalx = xbottom + 0.5*xheight*(1.0+
       cos(2.0*3.1415927/(xmax-xmin)*(x-(xmin+xmax)/2.0)));
    if (z <= zmin || z >= zmax) ylocalz = zbottom;
    else ylocalz = zbottom + 0.5*zheight*(1.0+
       cos(2.0*3.1415927/(zmax-zmin)*(z-(zmin+zmax)/2.0)));
    ylocal = ylocalx;
    if (ylocal > ylocalz) ylocal = ylocalz;
    if (y >= ylocal) isite = VACANT;
    else if (y > yfix) isite = OCCUPIED;
    else isite = FIXED;
    lattice[i] = isite;
  }
}

/* ---------------------------------------------------------------------- */

AppSurface::~AppSurface()
{
  delete random;
  delete [] sites;
  delete [] sitesSchwoebel;
  delete [] check;
  delete [] checkSchwoebel;
  memory->sfree(events);
  memory->sfree(firstevent);
  delete [] ecoord;
}

/* ---------------------------------------------------------------------- */

void AppSurface::init_app()
{
  // check used to temporarily label owned and ghost sites

  delete [] check;
  delete [] checkSchwoebel;
  check = new int[nlocal+nghost];
  checkSchwoebel = new int[nlocal+nghost];
  for (int i = 0; i < nlocal+nghost; i++) check[i] = checkSchwoebel[i] = 0;

  memory->sfree(events);
  memory->sfree(firstevent);

  events = NULL;
  nevents = maxevent = 0;
  firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;

  t_inverse = t_inverse/boltz;
}

/* ---------------------------------------------------------------------- */

void AppSurface::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"ecoord") == 0) {
    if (narg != 2) error->all("Illegal ecoord command");
    int index = atoi(arg[0]);
    double value = atof(arg[1]);
    if (index < 0 || index > maxneigh) error->all("Illegal ecoord command");
    ecoord[index] = value;
  } else if (strcmp(command,"vibrafreq") == 0) {
    if (narg != 1) error->all("Illegal vibrafreq command");
    vibrafreq = atof(arg[0]);
  } else if (strcmp(command,"ebarrier") == 0) {
    if (narg != 2) error->all("Illegal ebarrier command");
    ebarrier = atof(arg[0]);
    nminRegular = atoi(arg[1]);
  } else if (strcmp(command,"eSchwoebel") == 0) {
    if (narg != 3) error->all("Illegal eSchwoebel command");
    eSchwoebel = atof(arg[0]);
    nmaxSchwoebel = atoi(arg[1]);
    nminSchwoebel = atoi(arg[2]);
  } else if (strcmp(command,"bondener") == 0) {
    if (narg != 1) error->all("Illegal bondener command");
    bondener = atof(arg[0]);
  } else error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppSurface::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  if (isite != VACANT)  {
     for (int j = 0; j < numneigh[i]; j++)
        if (lattice[neighbor[i][j]] != VACANT) eng++;
  }
  return (double) 0.5*eng*bondener;

  /*
  int n = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) n++;
  return ecoord[n];
  */
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
------------------------------------------------------------------------- */

void AppSurface::site_event_rejection(int i, RandomPark *random)
{
  // event = exchange with random neighbor

  int iran = (int) (numneigh[i]*random->uniform());
  if (iran >= numneigh[i]) iran = numneigh[i] - 1;
  int j = neighbor[i][iran];

  double einitial = 2.0*site_energy(i);

  int mystate = lattice[i];
  int neighstate = lattice[j];
  lattice[i] = neighstate;
  lattice[j] = mystate;

  double efinal = 2.0*site_energy(j);

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

double AppSurface::site_propensity(int i)
{
  int ineigh,jneigh,j,k,nsites;

  //    Possible events:
  // 1: OCCUPIED site exchanges with an adjacent VACANT site.
  //    Do not allow the event if this vacant site has less than
  //    nminRegular OCCUPIED neighbors
  // 2: Schwoebel events = the occupant at site i jumps through a nearest
  //    VACANT site j1 to a second-nearest VACANT site k around a nearest
  //    OCCUPIED site j2. If the occupant initially has more than 
  //    nmaxSchwoebel nearest OCCUPIED neighbors or finally has less
  //    than nminSchwoebel OCCUPIED neighbors, do not allow this event.

  clear_events(i);

  int mystate = lattice[i];
  if (mystate == VACANT || mystate == FIXED) return 0.0;

  int neighstate,nbs;
  double proball = 0.0;
  double einitial,efinal,probone;

  einitial = 2.0*site_energy(i);
  for (ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    neighstate = lattice[j];
    if (neighstate == VACANT) {
      lattice[i] = neighstate;
      lattice[j] = mystate;
      efinal = 2.0*site_energy(j);
      nbs = (int) (efinal/bondener + 0.5);
      if (nbs >= nminRegular) {
        if (efinal <= einitial)
          probone = vibrafreq*exp(-ebarrier*t_inverse);
        else {
          probone = vibrafreq*
            exp((einitial-efinal-ebarrier)*t_inverse);
        }
        if (probone > 0.0) {
          add_event(i,j,probone,0);
          proball += probone;
        }
      }
      lattice[i] = mystate;
      lattice[j] = neighstate;
    }
  }

  if (eSchwoebel < 1.0e6) {

// flag the neighbors of i, VACANT = 1, OCCUPIED = 2
    for (ineigh = 0; ineigh < numneigh[i]; ineigh++) {
      j = neighbor[i][ineigh];
      if (lattice[j] == VACANT) {
        checkSchwoebel[j] = 1;
      } else {
        checkSchwoebel[j] = 2;
      }
    } 

    nbs = (int) (einitial/bondener + 0.5);
    if (nbs <= nmaxSchwoebel) {
// checkSchwoebel_k = 1 or 2 means i's neighbor, skip
// checkSchwoebel_k = 30 means having j1 and j2 as neighbors,
//                    skip because already seen
// checkSchwoebel_k = 10*checkSchwoebel_j means having detected
//                    j1 or j2 neighbor, but not yet detected
//                    another j1 or j2, skip
      nsites = 0;
      for (ineigh = 0; ineigh < numneigh[i]; ineigh++) {
        j = neighbor[i][ineigh];
        for (jneigh = 0; jneigh < numneigh[j]; jneigh++) {
          k = neighbor[j][jneigh];
          if (lattice[k] != VACANT) continue;
          if (checkSchwoebel[k] == 1) continue;
          if (checkSchwoebel[k] == 2) continue;
          if (checkSchwoebel[k] == 30) continue;
          if (checkSchwoebel[k] == 10*checkSchwoebel[j]) continue;
          if (checkSchwoebel[k] == 0) sitesSchwoebel[nsites++] = k;
          checkSchwoebel[k] = checkSchwoebel[k] + 10*checkSchwoebel[j];
          if (checkSchwoebel[k] == 30) {
            neighstate = lattice[k];
            lattice[i] = neighstate;
            lattice[k] = mystate;
            efinal = 2.0*site_energy(k);
            nbs = (int) (efinal/bondener + 0.5);
            if (nbs >= nminSchwoebel) {
              if (efinal <= einitial)
                probone = vibrafreq*exp(-eSchwoebel*t_inverse);
              else {
                probone = vibrafreq*
                  exp((einitial-efinal-eSchwoebel)*t_inverse);
               }
              if (probone > 0.0) {
                add_event(i,k,probone,1);
                proball += probone;
              }
            }
            lattice[i] = mystate;
            lattice[k] = neighstate;
          }
        }
      }
      for (j = 0; j < nsites; j++) {
        k = sitesSchwoebel[j];
        checkSchwoebel[k] = 0;
      }
    }
    for (ineigh = 0; ineigh < numneigh[i]; ineigh++) {
      j = neighbor[i][ineigh];
      checkSchwoebel[j] = 0;
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   ignore neighbor sites that should not be updated (isite < 0)
------------------------------------------------------------------------- */

void AppSurface::site_event(int i, class RandomPark *random)
{
  int j,jj,k,kk,kkk,m,mm,mmm,isite,jumptype;

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
  jumptype = events[ievent].Schwoebel;
  lattice[i] = VACANT;
  lattice[j] = OCCUPIED;

  // compute propensity changes for self and swap site and their 1,2,3 neighs
  // depending on Schwoebel barrier
  // use check[] to avoid resetting propensity of same site

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  int nsites = 0;
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
    if (isite < 0) continue;
    if (check[isite] == 0) {
      sites[nsites++] = isite;
      propensity[isite] = site_propensity(m);
      check[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite < 0) continue;
      if (check[isite] == 0) {
        sites[nsites++] = isite;
        propensity[isite] = site_propensity(mm);
        check[isite] = 1;
      }
//      if (jumptype == 1) {
      if (eSchwoebel < 1.0e6) {
        for (kkk = 0; kkk < numneigh[mm]; kkk++) {
          mmm = neighbor[mm][kkk];
          isite = i2site[mmm];
          if (isite < 0) continue;
          if (check[isite] == 0) {
            sites[nsites++] = isite;
            propensity[isite] = site_propensity(mmm);
            check[isite] = 1;
          }
        }
      }
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
      if (isite < 0) continue;
      if (check[isite] == 0) {
        sites[nsites++] = isite;
        propensity[isite] = site_propensity(mm);
        check[isite] = 1;
      }
//      if (jumptype == 1) {
      if (eSchwoebel < 1.0e6) {
        for (kkk = 0; kkk < numneigh[mm]; kkk++) {
          mmm = neighbor[mm][kkk];
          isite = i2site[mmm];
          if (isite < 0) continue;
          if (check[isite] == 0) {
            sites[nsites++] = isite;
            propensity[isite] = site_propensity(mmm);
            check[isite] = 1;
          }
        }
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

void AppSurface::clear_events(int i)
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

void AppSurface::add_event(int i, int partner, double propensity, int jumptype)
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
  events[freeevent].Schwoebel = jumptype;
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}
