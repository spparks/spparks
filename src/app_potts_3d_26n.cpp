/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_3d_26n.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPotts3d26n::AppPotts3d26n(SPPARKS *spk, int narg, char **arg) : 
  AppPotts3d(spk,narg,arg)
{
  // parse any remaining arguments

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sample_argument") == 0) {
      iarg ++;
    } else {
      error->all("Illegal app_style potts/3d/26n command");
    }
  }

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_3d_T_array(lattice,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
			    "app:lattice");

  // initialize my portion of lattice

  if (init_style == RANDOM) {
  // each site = one of nspins
  // loop over global list so assignment is independent of # of procs

    int i,j,k,ii,jj,kk,isite;
    for (i = 1; i <= nx_global; i++) {
      ii = i - nx_offset;
      for (j = 1; j <= ny_global; j++) {
	jj = j - ny_offset;
	for (k = 1; k <= nz_global; k++) {
	  kk = k - nz_offset;
	  isite = random->irandom(nspins);
	  if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local &&
	      kk >= 1 && kk <= nz_local) {
	    lattice[ii][jj][kk] = isite;
	  }
	}
      } 
    }
  } else if (init_style == READ) read_spins(spinfile);
}

/* ---------------------------------------------------------------------- */

AppPotts3d26n::~AppPotts3d26n()
{
  memory->destroy_3d_T_array(lattice);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts3d26n::site_energy(int i, int j, int k)
{
  int ii,jj,kk;

  int isite = lattice[i][j][k];
  int eng = 0;
  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	if (isite != lattice[ii][jj][kk]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
   if site changes, unset mask of all neighbor sites with affected propensity
------------------------------------------------------------------------- */

void AppPotts3d26n::site_event_rejection(int i, int j, int k,
					RandomPark *random)
{
  int oldstate = lattice[i][j][k];
  double einitial = site_energy(i,j,k);

  // event = random spin

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  lattice[i][j][k] = iran;
  double efinal = site_energy(i,j,k);

  // event = random neighbor spin

  //int iran = (int) (26.0*random->uniform()) + 1;
  //if (iran > 26) iran = 26;
  //if (iran > 13) iran++;
  //int ii = iran % 3;
  //int jj = (iran/3) % 3;
  //int kk = iran / 9;
  //return lattice[i-1+ii][j-1+jj][k-1+kk];

  // event = random unique neighbor spin
  // not yet implemented

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i][j][k] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i][j][k] = oldstate;
  }

  if (Lmask) {
    if (einitial < 3.0) mask[i][j][k] = 1;
    if (lattice[i][j][k] != oldstate) {
      int ii,jj,kk;
      for (ii = i-1; ii <= i+1; ii++)
	for (jj = j-1; jj <= j+1; jj++)
	  for (kk = k-1; kk <= k+1; kk++)
	    mask[ii][jj][kk] = 0;
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

double AppPotts3d26n::site_propensity(int i, int j, int k)
{
  // possible events = spin flips to neighboring site different than self

  int sites[27];
  int oldstate = lattice[i][j][k];
  int nevent = 0;

  int ii,jj,kk;
  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	add_unique(oldstate,nevent,sites,ii,jj,kk);

  // for each possible flip:
  // compute energy difference between initial and final state
  // sum to prob for all events on this site

  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < nevent; m++) {
    lattice[i][j][k] = sites[m];
    efinal = site_energy(i,j,k);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
  }

  lattice[i][j][k] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, adjust neighbor sites indices for PBC
   if proc owns sector, ignore non-updated neighbors (isite < 0)
------------------------------------------------------------------------- */

void AppPotts3d26n::site_event(int i, int j, int k,
			      int full, RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[ijk2site[i][j][k]];

  // possible events = spin flips to neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int sites[27];
  int oldstate = lattice[i][j][k];
  int nevent = 0;

  int ii,jj,kk;
  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	add_unique(oldstate,nevent,sites,ii,jj,kk);

  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < nevent; m++) {
    lattice[i][j][k] = sites[m];
    efinal = site_energy(i,j,k);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  if (full) update_ghost_sites(i,j,k);

  // compute propensity changes for self and neighbor sites

  int iloop,jloop,kloop,isite;

  int nsites = 0;

  for (iloop = i-2; iloop <= i+2; iloop++)
    for (jloop = j-2; jloop <= j+2; jloop++)
      for (kloop = k-2; kloop <= k+2; kloop++) {
	ii = iloop; jj = jloop; kk = kloop;
	if (full) ijkpbc(ii,jj,kk);
	isite = ijk2site[ii][jj][kk];
	if (isite >= 0) {
	  sites[nsites++] = isite;
	  propensity[isite] = site_propensity(ii,jj,kk);
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

void AppPotts3d26n::push_connected_neighbors(int i, int j, int k, 
					     int*** cluster_ids, int id, 
					     std::stack<int>* cluststack)
{
  int iii,jjj,kkk;

  for (int ii = i-1; ii <= i+1; ii++) {
    for (int jj = j-1; jj <= j+1; jj++) {
      for (int kk = k-1; kk <= k+1; kk++) {
	if (lattice[ii][jj][kk] == lattice[i][j][k] &&
	    cluster_ids[ii][jj][kk] == 0) {
	  cluststack->push(ii);
	  cluststack->push(jj);
	  cluststack->push(kk);
	  cluster_ids[ii][jj][kk] = id;
	}
      }
    }
  }

}

/* ----------------------------------------------------------------------
   Add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */

void AppPotts3d26n::connected_ghosts(int i, int j, int k, 
				     int*** cluster_ids, Cluster* clustlist, 
				     int idoffset)
{
  int iclust;

  for (int ii = i-1; ii <= i+1; ii++) {
    for (int jj = j-1; jj <= j+1; jj++) {
      for (int kk = k-1; kk <= k+1; kk++) {
	if (lattice[ii][jj][kk] == lattice[i][j][k] &&
	    (ii < 1 || ii > nx_local || 
	     jj < 1 || jj > ny_local ||
	     kk < 1 || kk > nz_local )) {
	  iclust = cluster_ids[i][j][k]-idoffset;
	  // Add ghost cluster to neighbors of local cluster
	  clustlist[iclust].add_neigh(cluster_ids[ii][jj][kk]);
	}
      }
    }
  }

}

