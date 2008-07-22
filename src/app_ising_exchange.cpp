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
#include "app_ising_exchange.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppIsingExchange::AppIsingExchange(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;

  // parse arguments

  if (narg < 2) error->all("Illegal app_style command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);

  options(narg-2,&arg[2]);

  // define lattice and partition it across processors

  create_lattice();
  sites = new int[1 + maxneigh];

  // initialize my portion of lattice
  // each site = one of 2 spins
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  std::map<int,int> hash;
  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));
  std::map<int,int>::iterator loc;

  int ilocal,isite;
  for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
    isite = random->irandom(2);
    loc = hash.find(iglobal);
    if (loc != hash.end()) lattice[loc->second] = isite;
  }
}

/* ---------------------------------------------------------------------- */

AppIsingExchange::~AppIsingExchange()
{
  delete random;
  delete [] sites;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsingExchange::site_energy(int i)
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

void AppIsingExchange::site_event_rejection(int i, RandomPark *random)
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

double AppIsingExchange::site_propensity(int i)
{
  // possible events = exchange with neighboring site different than self

  int mystate = lattice[i];

  int neighstate;
  double einitial,efinal;
  double prob = 0.0;

  for (int j = 0; j < numneigh[i]; j++) {
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

void AppIsingExchange::site_event(int i, class RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[i2site[i]];

  // possible events = exchange with neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int mystate = lattice[i];

  int neighstate;
  double einitial,efinal;
  double prob = 0.0;

  for (int j = 0; j < numneigh[i]; j++) {
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

  // compute propensity changes for self and neighbor sites

  int m;

  int nsites = 0;
  int isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);

  for (int j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
  }

  solve->update(nsites,sites,propensity);
}
