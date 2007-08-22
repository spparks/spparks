/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppPotts::AppPotts(SPK *spk, int narg, char **arg) : AppLattice(spk,narg,arg)
{
  // parse arguments

  if (narg < 3) error->all("Invalid app_style potts command");

  seed = atoi(arg[1]);
  nspins = atoi(arg[2]);

  random = new RandomPark(seed);

  options(narg-3,&arg[3]);

  // define lattice and partition it across processors
  
  create_lattice();
  lattice = (int *) memory->smalloc((nlocal+nghost)*sizeof(int),"app:lattice");
  sites = new int[1 + maxconnect];

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assigment is independent of # of procs

  int ilocal,isite;
  for (int iglobal = 0; iglobal < nglobal; iglobal++) {
    isite = random->irandom(nspins);
    ilocal = iglobal;     // change this line in parallel to a map lookup
    lattice[ilocal] = isite;
  }

  // setup communicator for ghost sites

  comm = new CommLattice(spk);
  if (dimension == 2)
    comm->init(nlocal,procwest,proceast,procsouth,procnorth,
	       delghost,dellocal);
  else if (dimension == 3)
    comm->init(nlocal,procwest,proceast,procsouth,procnorth,procdown,procup,
	       delghost,dellocal);
}

/* ---------------------------------------------------------------------- */

AppPotts::~AppPotts()
{
  delete random;
  memory->sfree(lattice);
  delete [] sites;
  delete comm;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != lattice[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site
------------------------------------------------------------------------- */

int AppPotts::site_pick_random(int i, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

int AppPotts::site_pick_local(int i, double ran)
{
  int iran = (int) (numneigh[i]*ran) + 1;
  if (iran > numneigh[i]) iran = numneigh[i];
  return lattice[neighbor[i][iran]];
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   full flag is ignored, since if 1 proc owns full domain, there are no ghosts
------------------------------------------------------------------------- */

double AppPotts::site_propensity(int i, int full)
{
  int j,k,value;
  double efinal;

  // possible events = spin flips to neighboring site different than self

  int nevent = 0;
  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (k = 0; k < nevent; k++)
      if (value == sites[k]) break;
    if (k < nevent) continue;
    sites[nevent++] = value;
  }

  // for each possible flip:
  // compute energy difference between initial and final state
  // sum to prob for all events on this site

  int oldstate = lattice[i];
  double einitial = site_energy(i);
  double prob = 0.0;

  for (k = 0; k < nevent; k++) {
    lattice[i] = sites[k];
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
  }

  lattice[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   full flag is ignored, since if 1 proc owns full domain, there are no ghosts
   don't compute propensities of ghost sites
------------------------------------------------------------------------- */

void AppPotts::site_event(int i, int full)
{
  int j,k,isite,value;
  double efinal;

  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[i];

  // possible events = spin flips to neighboring site different than self
  // find one event, accumulate its probability
  // compare prob to threshhold, break when reach it to select event

  double einitial = site_energy(i);
  double prob = 0.0;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (k = 0; k < nevent; k++)
      if (value == sites[k]) break;
    if (k < nevent) continue;
    sites[nevent++] = value;

    lattice[i] = sites[k];
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);

    if (prob >= threshhold) break;
  }

  // compute propensity changes for self and neighbor sites

  int nsites = 0;
  sites[nsites++] = i;
  propensity[i] = site_propensity(i,full);

  for (j = 0; j < numneigh[i]; j++) {
    isite = neighbor[i][j];
    if (isite >= nlocal) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(isite,full);
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
  OK to clear ghost site mask values
------------------------------------------------------------------------- */

void AppPotts::site_clear_mask(char *mask, int i)
{
  mask[i] = 0;
  for (int j = 0; j < numneigh[i]; j++) mask[neighbor[i][j]] = 0;
}
