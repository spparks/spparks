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
#include "stdlib.h"
#include "string.h"
#include "app_erbium.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NOOP,FCC,OCTA,TETRA};
enum{ZERO,ERBIUM,HYDROGEN,HELIUM,VACANCY};      // same as DiagErbium

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppErbium::AppErbium(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = trate = NULL;
  spropensity = dpropensity = tpropensity = NULL;
  stype = sinput = soutput = NULL;
  dtype = dinput = doutput = NULL;
  ttype = tinput = toutput = NULL;
  scount = dcount = tcount = NULL;
}

/* ---------------------------------------------------------------------- */

AppErbium::~AppErbium()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);

  memory->destroy(srate);
  memory->destroy(drate);
  memory->destroy(trate);
  memory->destroy(spropensity);
  memory->destroy(dpropensity);
  memory->destroy(tpropensity);
  memory->destroy(stype);
  memory->destroy(sinput);
  memory->destroy(soutput);
  memory->destroy(dtype);
  memory->destroy(dinput);
  memory->destroy(doutput);
  memory->destroy(ttype);
  memory->destroy(tinput);
  memory->destroy(toutput);
  memory->destroy(scount);
  memory->destroy(dcount);
  memory->destroy(tcount);
}

/* ---------------------------------------------------------------------- */

void AppErbium::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 5) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) stype[none] = FCC;
      else if (strcmp(arg[1],"oct") == 0) stype[none] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) stype[none] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"er") == 0) sinput[none] = ERBIUM;
      else if (strcmp(arg[2],"h") == 0) sinput[none] = HYDROGEN;
      else if (strcmp(arg[2],"he") == 0) sinput[none] = HELIUM;
      else if (strcmp(arg[2],"vac") == 0) sinput[none] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      srate[none] = atof(arg[3]);

      if (strcmp(arg[4],"er") == 0) soutput[none] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) soutput[none] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) soutput[none] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) soutput[none] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      none++;
      
    } else if (rstyle == 2) {
      if (narg != 8) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) dtype[ntwo][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) dtype[ntwo][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) dtype[ntwo][0] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) dtype[ntwo][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) dtype[ntwo][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) dtype[ntwo][1] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"er") == 0) dinput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[3],"h") == 0) dinput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[3],"he") == 0) dinput[ntwo][0] = HELIUM;
      else if (strcmp(arg[3],"vac") == 0) dinput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"er") == 0) dinput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) dinput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) dinput[ntwo][1] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) dinput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      drate[ntwo] = atof(arg[5]);

      if (strcmp(arg[6],"er") == 0) doutput[ntwo][0] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) doutput[ntwo][0] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) doutput[ntwo][0] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) doutput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"er") == 0) doutput[ntwo][1] = ERBIUM;
      else if (strcmp(arg[7],"h") == 0) doutput[ntwo][1] = HYDROGEN;
      else if (strcmp(arg[7],"he") == 0) doutput[ntwo][1] = HELIUM;
      else if (strcmp(arg[7],"vac") == 0) doutput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      ntwo++;

    } else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"fcc") == 0) ttype[nthree][0] = FCC;
      else if (strcmp(arg[1],"oct") == 0) ttype[nthree][0] = OCTA;
      else if (strcmp(arg[1],"tet") == 0) ttype[nthree][0] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"fcc") == 0) ttype[nthree][1] = FCC;
      else if (strcmp(arg[2],"oct") == 0) ttype[nthree][1] = OCTA;
      else if (strcmp(arg[2],"tet") == 0) ttype[nthree][1] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"fcc") == 0) ttype[nthree][2] = FCC;
      else if (strcmp(arg[3],"oct") == 0) ttype[nthree][2] = OCTA;
      else if (strcmp(arg[3],"tet") == 0) ttype[nthree][2] = TETRA;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"er") == 0) tinput[nthree][0] = ERBIUM;
      else if (strcmp(arg[4],"h") == 0) tinput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[4],"he") == 0) tinput[nthree][0] = HELIUM;
      else if (strcmp(arg[4],"vac") == 0) tinput[nthree][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[5],"er") == 0) tinput[nthree][1] = ERBIUM;
      else if (strcmp(arg[5],"h") == 0) tinput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[5],"he") == 0) tinput[nthree][1] = HELIUM;
      else if (strcmp(arg[5],"vac") == 0) tinput[nthree][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[6],"er") == 0) tinput[nthree][2] = ERBIUM;
      else if (strcmp(arg[6],"h") == 0) tinput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[6],"he") == 0) tinput[nthree][2] = HELIUM;
      else if (strcmp(arg[6],"vac") == 0) tinput[nthree][2] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      trate[nthree] = atof(arg[7]);

      if (strcmp(arg[8],"er") == 0) toutput[nthree][0] = ERBIUM;
      else if (strcmp(arg[8],"h") == 0) toutput[nthree][0] = HYDROGEN;
      else if (strcmp(arg[8],"he") == 0) toutput[nthree][0] = HELIUM;
      else if (strcmp(arg[8],"vac") == 0) toutput[nthree][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[9],"er") == 0) toutput[nthree][1] = ERBIUM;
      else if (strcmp(arg[9],"h") == 0) toutput[nthree][1] = HYDROGEN;
      else if (strcmp(arg[9],"he") == 0) toutput[nthree][1] = HELIUM;
      else if (strcmp(arg[9],"vac") == 0) toutput[nthree][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[10],"er") == 0) toutput[nthree][2] = ERBIUM;
      else if (strcmp(arg[10],"h") == 0) toutput[nthree][2] = HYDROGEN;
      else if (strcmp(arg[10],"he") == 0) toutput[nthree][2] = HELIUM;
      else if (strcmp(arg[10],"vac") == 0) toutput[nthree][2] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      nthree++;

    } else error->all(FLERR,"Illegal event command");
  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppErbium::grow_app()
{
  type = iarray[0];
  element = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppErbium::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = new int[3 + 3*maxneigh];
  }

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] < FCC || type[i] > TETRA) flag = 1;
    if (element[i] < ERBIUM || element[i] > VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- */

void AppErbium::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates

  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app erbium");

  for (int m = 0; m < none; m++) {
    spropensity[m] = srate[m];
    scount[m] = 0;
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = exp(-drate[m]/temperature);
    dcount[m] = 0;
  }
  for (int m = 0; m < nthree; m++) {
    tpropensity[m] = exp(-trate[m]/temperature);
    tcount[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppErbium::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppErbium::site_propensity(int i)
{
  int j,k,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;

  // currently fcc sites have no events
  // can remove this later

  if (type[i] == FCC) return proball = 0.0;

  // single-site events

  for (m = 0; m < none; m++) {
    if (type[i] != stype[m] || element[i] != sinput[m]) continue;
    add_event(i,1,m,spropensity[m],-1,-1);
    proball += spropensity[m];
  }

  // double-site events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ntwo; m++) {
      if (type[i] != dtype[m][0] || element[i] != dinput[m][0]) continue;
      if (type[j] != dtype[m][1] || element[j] != dinput[m][1]) continue;
      add_event(i,2,m,dpropensity[m],j,-1);
      proball += dpropensity[m];
    }
  }

  // triple-site events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (int kk = 0; kk < numneigh[i]; kk++) {
      if (jj == kk) continue;
      k = neighbor[i][kk];
      for (m = 0; m < nthree; m++) {
	if (type[i] != ttype[m][0] || element[i] != tinput[m][0]) continue;
	if (type[j] != ttype[m][1] || element[j] != tinput[m][1]) continue;
	if (type[k] != ttype[m][2] || element[k] != tinput[m][2]) continue;
	add_event(i,3,m,tpropensity[m],j,k);
	proball += tpropensity[m];
      }
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppErbium::site_event(int i, class RandomPark *random)
{
  int j,k,m,n;

  // pick one event from total propensity by accumulating its probability
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

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;

  int ibef = 0;
  for (int jj = 0; jj < nlocal; jj++)
    if (element[jj] == HYDROGEN) ibef++;

  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
  } else if (rstyle == 2) {
    element[i] = doutput[which][0];
    element[j] = doutput[which][1];
    dcount[which]++;
  } else {
    element[i] = toutput[which][0];
    element[j] = toutput[which][1];
    element[k] = toutput[which][2];
    tcount[which]++;
  }

  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (rstyle >= 2) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(m);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  if (rstyle == 3) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(m);
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

void AppErbium::add_event(int i, int rstyle, int which, double propensity,
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

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single, double, or triple
------------------------------------------------------------------------- */

void AppErbium::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    memory->grow(srate,n,"app/erbium:srate");
    memory->grow(spropensity,n,"app/erbium:spropensity");
    memory->grow(stype,n,"app/erbium:stype");
    memory->grow(sinput,n,"app/erbium:sinput");
    memory->grow(soutput,n,"app/erbium:soutput");
    memory->grow(scount,n,"app/erbium:scount");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    memory->grow(drate,n,"app/erbium:drate");
    memory->grow(dpropensity,n,"app/erbium:dpropensity");
    dtype = memory->grow(dtype,n,2,"app/erbium:dtype");
    dinput = memory->grow(dinput,n,2,"app/erbium:dinput");
    doutput = memory->grow(doutput,n,2,"app/erbium:doutput");
    memory->grow(dcount,n,"app/erbium:dcount");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    memory->grow(trate,n,"app/erbium:trate");
    memory->grow(tpropensity,n,"app/erbium:tpropensity");
    ttype = memory->grow(ttype,n,3,"app/erbium:ttype");
    tinput = memory->grow(tinput,n,3,"app/erbium:tinput");
    toutput = memory->grow(toutput,n,3,"app/erbium:toutput");
    memory->grow(tcount,n,"app/erbium:tcount");
  }
}
