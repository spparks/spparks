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
#include "app_potts_pin.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsPin::AppPottsPin(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg) {}

/* ---------------------------------------------------------------------- */

void AppPottsPin::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"pin") == 0) {
    if (narg != 1) error->all("Illegal pin command");
    pfraction = atof(arg[0]);
    if (pfraction < 0.0 || pfraction > 1.0)
      error->all("Illegal pin command");
    pin_create();
  } else error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins, but not to a pinned site
------------------------------------------------------------------------- */

void AppPottsPin::site_event_rejection(int i, RandomPark *random)
{
  // no events for a pinned site

  if (lattice[i] > nspins) return;

  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // event = random spin from 1 to nspins, including self

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  lattice[i] = iran;
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion
  // null bin extends to nspins

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  if (lattice[i] != oldstate) naccept++;

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (lattice[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppPottsPin::site_propensity(int i)
{
  // no events for a pinned site

  if (lattice[i] > nspins) return 0.0;

  // events = spin flips to neighboring site different than self
  // disallow flip to pinned site
  // disallow wild flips = flips to value different than all neighs

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i] || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  // for each flip:
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  int oldstate = lattice[i];
  double einitial = site_energy(i);
  double efinal;
  double prob = 0.0;

  for (m = 0; m < nevent; m++) {
    lattice[i] = unique[m];
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
  }

  lattice[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPottsPin::site_event(int i, RandomPark *random)
{
  int j,m,value;

  // pick one event from total propensity by accumulating its probability
  // disallow flip to pinned site
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double efinal;

  int oldstate = lattice[i];
  double einitial = site_energy(i);
  double prob = 0.0;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == oldstate || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;

    lattice[i] = value;
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  // compute propensity changes for self and neighbor sites
  // ignore update of neighbor sites with isite < 0

  int nsites = 0;
  int isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   change pfraction of sites to pinned sites
------------------------------------------------------------------------- */

void AppPottsPin::pin_create()
{
  int i,nattempt,iglobal,nme,npin;

  int ndesired = static_cast<int> (pfraction*nglobal);

  RandomPark *random = new RandomPark(ranmaster->uniform());

  std::map<int,int> hash;
  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));
  std::map<int,int>::iterator loc;

  npin = 0;
  while (npin < ndesired) {
    nattempt = ndesired - npin;
    for (i = 0; i < nattempt; i++) {
      iglobal = random->irandom(nglobal);
      loc = hash.find(iglobal);
      if (loc != hash.end()) lattice[loc->second] = nspins+1;
    }

    nme = 0;
    for (i = 0; i < nlocal; i++)
      if (lattice[i] > nspins) nme++;
    MPI_Allreduce(&nme,&npin,1,MPI_INT,MPI_SUM,world);
  }

  delete random;
}
