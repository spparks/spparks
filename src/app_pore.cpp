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
#include "app_pore.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPore::AppPore(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;

  // parse arguments

  if (narg < 7) error->all("Illegal app_style command");

  double xc = atof(arg[1]);
  double yc = atof(arg[2]);
  double zc = atof(arg[3]);
  double diameter = atof(arg[4]);
  double thickness = atof(arg[5]);
  int seed = atoi(arg[6]);
  random = new RandomPark(seed);

  options(narg-7,&arg[7]);

  // define lattice and partition it across processors
  // sites must be large enough for 2 sites and their 1st/2nd nearest neighbors

  create_lattice();
  sites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  check = NULL;

  // initialize my portion of lattice
  // each site = 1 (vacancy) or 2 (occupied)
  // pore geometry defines occupied vs unoccupied

  double x,y,z;
  int isite;
  for (int i = 0; i < nlocal; i++) {
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    if (z > zc + 0.5*thickness || z < zc - 0.5*thickness) isite = 1;
    else isite = 2;
    if (isite == 2) {
      if ((x-xc)*(x-xc) + (y-yc)*(y-yc) < 0.25*diameter*diameter) isite = 1;
    }
    lattice[i] = isite;
  }
}

/* ---------------------------------------------------------------------- */

AppPore::~AppPore()
{
  delete random;
  delete [] sites;
  delete [] check;
}

/* ---------------------------------------------------------------------- */

void AppPore::init_app()
{
  delete [] check;
  check = new int[nlocal];
  for (int i = 0; i < nlocal; i++) check[i] = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPore::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != lattice[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
------------------------------------------------------------------------- */

void AppPore::site_event_rejection(int i, RandomPark *random)
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
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, there are no ghosts, so ignore full flag
------------------------------------------------------------------------- */

double AppPore::site_propensity(int i)
{
  int j;

  // possible events = exchange with neighboring site different than self

  int mystate = lattice[i];

  int neighstate;
  double einitial,efinal;
  double prob = 0.0;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    neighstate = lattice[j];
    if (neighstate != mystate) {
      einitial = site_energy(i) + site_energy(j);
      lattice[i] = neighstate;
      lattice[j] = mystate;
      efinal = site_energy(i) + site_energy(j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
      lattice[i] = mystate;
      lattice[j] = neighstate;
    }
  }

  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, there are no ghosts, so ignore full flag
   if proc owns sector, ignore neighbor sites that are ghosts
------------------------------------------------------------------------- */

void AppPore::site_event(int i, class RandomPark *random)
{
  int j,jj,k,kk,m,mm;

  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[i2site[i]];

  // possible events = exchange with neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int mystate = lattice[i];

  int neighstate;
  double einitial,efinal;
  double prob = 0.0;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    neighstate = lattice[j];
    if (neighstate != mystate) {
      einitial = site_energy(i) + site_energy(j);
      lattice[i] = neighstate;
      lattice[j] = mystate;
      efinal = site_energy(i) + site_energy(j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
      if (prob >= threshhold) break;
      lattice[i] = mystate;
      lattice[j] = neighstate;
    }
  }

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
