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
#include "app_potts_2d_8n.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "cluster.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppPotts2d8n::AppPotts2d8n(SPPARKS *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Illegal app_style command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  int seed = atoi(arg[4]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nxlo,nxhi,nylo,nyhi,
			    "applattice2d:lattice");

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assignment is independent of # of procs

  int i,j,ii,jj,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      isite = random->irandom(nspins);
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local) {
	lattice[ii][jj] = isite;
      }
    } 
  }
}

/* ---------------------------------------------------------------------- */

AppPotts2d8n::~AppPotts2d8n()
{
  delete random;
  memory->destroy_2d_T_array(lattice,nxlo,nylo);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts2d8n::site_energy(int i, int j)
{
  int isite = lattice[i][j];
  int eng = 0;
  if (isite != lattice[i-1][j-1]) eng++;
  if (isite != lattice[i-1][j]) eng++;
  if (isite != lattice[i-1][j+1]) eng++;
  if (isite != lattice[i][j-1]) eng++;
  if (isite != lattice[i][j+1]) eng++;
  if (isite != lattice[i+1][j-1]) eng++;
  if (isite != lattice[i+1][j]) eng++;
  if (isite != lattice[i+1][j+1]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
   if site changes, unset mask of all neighbor sites with affected propensity
------------------------------------------------------------------------- */

void AppPotts2d8n::site_event_rejection(int i, int j, RandomPark *random)
{
  int oldstate = lattice[i][j];
  double einitial = site_energy(i,j);

  // event = random spin

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  lattice[i][j] = iran;
  double efinal = site_energy(i,j);

  // event = random neighbor spin

  //int iran = (int) (8.0*random->uniform());
  //if (iran == 0) return lattice[i-1][j-1];
  //else if (iran == 1) return lattice[i-1][j];
  //else if (iran == 2) return lattice[i-1][j+1];
  //else if (iran == 3) return lattice[i][j-1];
  //else if (iran == 4) return lattice[i][j+1];
  //else if (iran == 5) return lattice[i+1][j-1];
  //else if (iran == 6) return lattice[i+1][j];
  //else return lattice[i+1][j+1];

  // event = random unique neighbor spin
  // not yet implemented

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i][j] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i][j] = oldstate;
  }

  if (Lmask) {
    if (einitial < 4.0) mask[i][j] = 1;
    if (lattice[i][j] != oldstate) {
      mask[i-1][j-1] = 0;
      mask[i-1][j] = 0;
      mask[i-1][j+1] = 0;
      mask[i][j-1] = 0;
      mask[i][j+1] = 0;
      mask[i+1][j-1] = 0;
      mask[i+1][j] = 0;
      mask[i+1][j+1] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
   propensity for one event is based on einitial,efinal
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity = Boltzmann factor
------------------------------------------------------------------------- */

double AppPotts2d8n::site_propensity(int i, int j)
{
  // possible events = spin flips to neighboring site different than self

  int sites[9];
  int oldstate = lattice[i][j];
  int nevent = 0;

  add_unique(oldstate,nevent,sites,i-1,j-1);
  add_unique(oldstate,nevent,sites,i-1,j);
  add_unique(oldstate,nevent,sites,i-1,j+1);
  add_unique(oldstate,nevent,sites,i,j-1);
  add_unique(oldstate,nevent,sites,i,j+1);
  add_unique(oldstate,nevent,sites,i+1,j-1);
  add_unique(oldstate,nevent,sites,i+1,j);
  add_unique(oldstate,nevent,sites,i+1,j+1);

  // for each possible flip:
  // compute energy difference between initial and final state
  // sum to prob for all events on this site

  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < nevent; m++) {
    lattice[i][j] = sites[m];
    efinal = site_energy(i,j);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
  }

  lattice[i][j] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, adjust neighbor sites indices for PBC
   if proc owns sector, ignore non-updated neighbors (isite < 0)
------------------------------------------------------------------------- */

void AppPotts2d8n::site_event(int i, int j, int full, RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[ij2site[i][j]];

  // possible events = spin flips to neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int sites[9];
  int oldstate = lattice[i][j];
  int nevent = 0;

  add_unique(oldstate,nevent,sites,i-1,j-1);
  add_unique(oldstate,nevent,sites,i-1,j);
  add_unique(oldstate,nevent,sites,i-1,j+1);
  add_unique(oldstate,nevent,sites,i,j-1);
  add_unique(oldstate,nevent,sites,i,j+1);
  add_unique(oldstate,nevent,sites,i+1,j-1);
  add_unique(oldstate,nevent,sites,i+1,j);
  add_unique(oldstate,nevent,sites,i+1,j+1);

  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < nevent; m++) {
    lattice[i][j] = sites[m];
    efinal = site_energy(i,j);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  if (full) update_ghost_sites(i,j);

  // compute propensity changes for self and neighbor sites

  int iloop,jloop,ii,jj,isite;

  int nsites = 0;

  for (iloop = i-1; iloop <= i+1; iloop++)
    for (jloop = j-1; jloop <= j+1; jloop++) {
      ii = iloop; jj = jloop;
      if (full) ijpbc(ii,jj);
      isite = ij2site[ii][jj];
      if (isite >= 0) {
	sites[nsites++] = isite;
	propensity[isite] = site_propensity(ii,jj);
      }
    }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */

void AppPotts2d8n::push_connected_neighbors(int i, int j, 
					    int **cluster_ids, int id,
					    std::stack<int>* cluststack)
{
  int iii,jjj;

  for (int ii = i-1; ii <= i+1; ii++) {
    for (int jj = j-1; jj <= j+1; jj++) {
      if (lattice[ii][jj] == lattice[i][j] &&
	  cluster_ids[ii][jj] == 0) {
	cluststack->push(ii);
	cluststack->push(jj);
	cluster_ids[ii][jj] = id;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */

void AppPotts2d8n::connected_ghosts(int i, int j, int **cluster_ids, 
				    Cluster* clustlist, int idoffset)
{
  int iclust;

  for (int ii = i-1; ii <= i+1; ii++) {
    for (int jj = j-1; jj <= j+1; jj++) {
      if (lattice[ii][jj] == lattice[i][j] &&
	  (ii < 1 || ii > nx_local || 
	   jj < 1 || jj > ny_local )) {
	iclust = cluster_ids[i][j]-idoffset;
	// Add ghost cluster to neighbors of local cluster
        clustlist[iclust].add_neigh(cluster_ids[ii][jj]);
      }
    }
  }
}

