/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_3d_12n.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPotts3d12n::AppPotts3d12n(SPPARKS *spk, int narg, char **arg) : 
  AppPotts3d(spk,narg,arg)
{
  delevent = 0;
  delpropensity = 2;

  // parse any remaining arguments

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sample_argument") == 0) {
      iarg ++;
    } else {
      error->all("Illegal app_style potts/3d/12n command");
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

AppPotts3d12n::~AppPotts3d12n()
{
  memory->destroy_3d_T_array(lattice,nxlo,nylo,nzlo);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts3d12n::site_energy(int i, int j, int k)
{
  int isite = lattice[i][j][k];
  int eng = 0;
  if (isite != lattice[i][j][k-2]) eng++;
  if (isite != lattice[i][j][k-1]) eng++;
  if (isite != lattice[i][j][k+1]) eng++;
  if (isite != lattice[i][j][k+2]) eng++;
  if (isite != lattice[i][j-2][k]) eng++;
  if (isite != lattice[i][j-1][k]) eng++;
  if (isite != lattice[i][j+1][k]) eng++;
  if (isite != lattice[i][j+2][k]) eng++;
  if (isite != lattice[i-2][j][k]) eng++;
  if (isite != lattice[i-1][j][k]) eng++;
  if (isite != lattice[i+1][j][k]) eng++;
  if (isite != lattice[i+2][j][k]) eng++;
  return (double) eng;
}


/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
   if site changes, unset mask of all neighbor sites with affected propensity
------------------------------------------------------------------------- */

void AppPotts3d12n::site_event_rejection(int i, int j, int k,
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

  //int iran = (int) (12.0*random->uniform());
  //if (iran == 0) return lattice[i-2][j][k];
  //else if (iran == 1) return lattice[i-1][j][k];
  //else if (iran == 2) return lattice[i+1][j][k];
  //else if (iran == 3) return lattice[i+2][j][k];
  //else if (iran == 4) return lattice[i][j-2][k];
  //else if (iran == 5) return lattice[i][j-1][k];
  //else if (iran == 6) return lattice[i][j+1][k];
  //else if (iran == 7) return lattice[i][j+2][k];
  //else if (iran == 8) return lattice[i][j][k-2];
  //else if (iran == 9) return lattice[i][j][k-1];
  //else if (iran == 10) return lattice[i][j][k+1];
  //else return lattice[i][j][k+2];

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
    if (einitial < 6.0) mask[i][j][k] = 1;
    if (lattice[i][j][k] != oldstate) {
      mask[i-2][j][k] = 0;
      mask[i-1][j][k] = 0;
      mask[i+1][j][k] = 0;
      mask[i+2][j][k] = 0;
      mask[i][j-2][k] = 0;
      mask[i][j-1][k] = 0;
      mask[i][j+1][k] = 0;
      mask[i][j+2][k] = 0;
      mask[i][j][k-2] = 0;
      mask[i][j][k-1] = 0;
      mask[i][j][k+1] = 0;
      mask[i][j][k+2] = 0;
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

double AppPotts3d12n::site_propensity(int i, int j, int k)
{
  // possible events = spin flips to neighboring site different than self

  int sites[13];
  int oldstate = lattice[i][j][k];
  int nevent = 0;

  add_unique(oldstate,nevent,sites,i-2,j,k);
  add_unique(oldstate,nevent,sites,i-1,j,k);
  add_unique(oldstate,nevent,sites,i+1,j,k);
  add_unique(oldstate,nevent,sites,i+2,j,k);

  add_unique(oldstate,nevent,sites,i,j-2,k);
  add_unique(oldstate,nevent,sites,i,j-1,k);
  add_unique(oldstate,nevent,sites,i,j+1,k);
  add_unique(oldstate,nevent,sites,i,j+2,k);

  add_unique(oldstate,nevent,sites,i,j,k-2);
  add_unique(oldstate,nevent,sites,i,j,k-1);
  add_unique(oldstate,nevent,sites,i,j,k+1);
  add_unique(oldstate,nevent,sites,i,j,k+2);

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

void AppPotts3d12n::site_event(int i, int j, int k,
			      int full, RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[ijk2site[i][j][k]];

  // possible events = spin flips to neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int sites[13];
  int oldstate = lattice[i][j][k];
  int nevent = 0;

  add_unique(oldstate,nevent,sites,i-2,j,k);
  add_unique(oldstate,nevent,sites,i-1,j,k);
  add_unique(oldstate,nevent,sites,i+1,j,k);
  add_unique(oldstate,nevent,sites,i+2,j,k);

  add_unique(oldstate,nevent,sites,i,j-2,k);
  add_unique(oldstate,nevent,sites,i,j-1,k);
  add_unique(oldstate,nevent,sites,i,j+1,k);
  add_unique(oldstate,nevent,sites,i,j+2,k);

  add_unique(oldstate,nevent,sites,i,j,k-2);
  add_unique(oldstate,nevent,sites,i,j,k-1);
  add_unique(oldstate,nevent,sites,i,j,k+1);
  add_unique(oldstate,nevent,sites,i,j,k+2);

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

  int ii,jj,kk,isite;

  int nsites = 0;

  ii = i-2; jj = j; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i-1; jj = j; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i+1; jj = j; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i+2; jj = j; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j-2; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j-1; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j+1; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j+2; kk = k;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j; kk = k-2;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j; kk = k-1;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j; kk = k+1;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  ii = i; jj = j; kk = k+2;
  if (full) ijkpbc(ii,jj,kk);
  isite = ijk2site[ii][jj][kk];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,kk);
  }

  solve->update(nsites,sites,propensity);
}
