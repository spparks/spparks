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
#include "app_ising_2d_4n_exchange.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppIsing2d4nExchange::
AppIsing2d4nExchange(SPPARKS *spk, int narg, char **arg) :
  AppLattice2d(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;

  // parse arguments

  if (narg != 4) error->all("Illegal app_style command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  int seed = atoi(arg[3]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nxlo,nxhi,nylo,nyhi,
			    "app:lattice");

  // initialize my portion of lattice
  // each site = one of 2 spins
  // loop over global list so assigment is independent of # of procs

  int i,j,ii,jj,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      isite = random->irandom(2);
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local)
	lattice[ii][jj] = isite;
    }
  }
}

/* ---------------------------------------------------------------------- */

AppIsing2d4nExchange::~AppIsing2d4nExchange()
{
  delete random;
  memory->destroy_2d_T_array(lattice,nxlo,nylo);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsing2d4nExchange::site_energy(int i, int j)
{
  int isite = lattice[i][j];
  int eng = 0;
  if (isite != lattice[i][j-1]) eng++;
  if (isite != lattice[i][j+1]) eng++;
  if (isite != lattice[i-1][j]) eng++;
  if (isite != lattice[i+1][j]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
   if site changes, unset mask of all neighbor sites with affected propensity
------------------------------------------------------------------------- */

void AppIsing2d4nExchange::site_event_rejection(int i, int j,
						RandomPark *random)
{
  // event = exchange with random neighbor

  int ii,jj;
  int iran = (int) (4.0*random->uniform());
  if (iran == 0) {
    ii = i-1; jj = j;
  } else if (iran == 1) {
    ii = i+1; jj = j;
  } else if (iran == 2) {
    ii = i; jj = j-1;
  } else {
    ii = i; jj = j+1;
  }

  double einitial = site_energy(i,j) + site_energy(ii,jj);

  int mystate = lattice[i][j];
  int neighstate = lattice[ii][jj];
  lattice[i][j] = neighstate;
  lattice[ii][jj] = mystate;

  double efinal = site_energy(i,j) + site_energy(ii,jj);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  if (Lmask) {
    if (einitial < 2.0) mask[i][j] = 1;
    if (lattice[i][j] != mystate) {
      mask[i][j] = 0;
      mask[i-1][j] = 0;
      mask[i+1][j] = 0;
      mask[i][j-1] = 0;
      mask[i][j+1] = 0;

      mask[ii][jj] = 0;
      mask[ii-1][jj] = 0;
      mask[ii+1][jj] = 0;
      mask[ii][jj-1] = 0;
      mask[ii][jj+1] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
   propensity for one event based on einitial,efinal
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity = Boltzmann factor
------------------------------------------------------------------------- */

double AppIsing2d4nExchange::site_propensity(int i, int j)
{
  // possible events = exchange with neighboring site different than self

  int mystate = lattice[i][j];

  int ii,jj,neighstate;
  double einitial,efinal;
  double prob = 0.0;

  ii = i-1; jj = j;
  neighstate = lattice[ii][jj];
  if (neighstate != mystate) {
    einitial = site_energy(i,j) + site_energy(ii,jj);
    lattice[i][j] = neighstate;
    lattice[ii][jj] = mystate;
    double efinal = site_energy(i,j) + site_energy(ii,jj);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  ii = i+1; jj = j;
  neighstate = lattice[ii][jj];
  if (neighstate != mystate) {
    einitial = site_energy(i,j) + site_energy(ii,jj);
    lattice[i][j] = neighstate;
    lattice[ii][jj] = mystate;
    double efinal = site_energy(i,j) + site_energy(ii,jj);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  ii = i; jj = j-1;
  neighstate = lattice[ii][jj];
  if (neighstate != mystate) {
    einitial = site_energy(i,j) + site_energy(ii,jj);
    lattice[i][j] = neighstate;
    lattice[ii][jj] = mystate;
    double efinal = site_energy(i,j) + site_energy(ii,jj);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  ii = i; jj = j+1;
  neighstate = lattice[ii][jj];
  if (neighstate != mystate) {
    einitial = site_energy(i,j) + site_energy(ii,jj);
    lattice[i][j] = neighstate;
    lattice[ii][jj] = mystate;
    double efinal = site_energy(i,j) + site_energy(ii,jj);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, adjust neighbor sites indices for PBC
   if proc owns sector, ignore non-updated neighbors (isite < 0)
------------------------------------------------------------------------- */

void AppIsing2d4nExchange::site_event(int i, int j, int full,
				      RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[ij2site[i][j]];

  // possible events = exchange with neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int mystate = lattice[i][j];

  int ii,jj,neighstate;
  double einitial,efinal;
  double prob = 0.0;

  for (int m = 0; m < 4; m++) {
    if (m == 0) {
      ii = i-1; jj = j;
    } else if (m == 1) {
      ii = i+1; jj = j;
    } else if (m == 2) {
      ii = i; jj = j-1;
    } else {
      ii = i; jj = j+1;
    }
    neighstate = lattice[ii][jj];
    if (neighstate == mystate) continue;

    einitial = site_energy(i,j) + site_energy(ii,jj);
    lattice[i][j] = neighstate;
    lattice[ii][jj] = mystate;
    double efinal = site_energy(i,j) + site_energy(ii,jj);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
    lattice[i][j] = mystate;
    lattice[ii][jj] = neighstate;
  }

  if (full) update_ghost_sites(i,j);

  // compute propensity changes for self and neighbor sites

  int m = ii;
  int n = jj;
  int isite,sites[10];

  int nsites = 0;

  ii = i; jj = j;
  isite = ij2site[ii][jj];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(ii,jj);

  ii = i-1; jj = j;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i+1; jj = j;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i; jj = j-1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i; jj = j+1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = m; jj = n;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = m-1; jj = n;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = m+1; jj = n;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = m; jj = n-1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = m; jj = n+1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  solve->update(nsites,sites,propensity);
}
