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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppPotts3d26n::AppPotts3d26n(SPK *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
  // parse arguments

  if (narg != 6) error->all("Invalid app_style potts/3d/26n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nz_global = atoi(arg[3]);
  nspins = atoi(arg[4]);
  seed = atoi(arg[5]);
  random = new RandomPark(seed);

  masklimit = 13.0;

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_3d_T_array(lattice,nx_local+2,ny_local+2,nz_local+2,
			    "applattice3d:lattice");

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assigment is independent of # of procs

  int i,j,k,ii,jj,kk,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      for (k = 1; k <= nz_global; k++) {
	kk = k - nz_offset;
	isite = random->irandom(nspins);
	if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local &&
	    kk >= 1 && kk <= nz_local)
	  lattice[ii][jj][kk] = isite;
      }
    }
  }

  // setup communicator for ghost sites

  comm = new CommLattice3d(spk);
  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup);
}

/* ---------------------------------------------------------------------- */

AppPotts3d26n::~AppPotts3d26n()
{
  delete random;
  memory->destroy_3d_T_array(lattice);
  delete comm;
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
   randomly pick new state for site
------------------------------------------------------------------------- */

int AppPotts3d26n::site_pick_random(int i, int j, int k, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

int AppPotts3d26n::site_pick_local(int i, int j, int k, double ran)
{
  int iran = (int) (26*ran) + 1;
  if (iran > 26) iran = 26;
  if (iran > 13) iran++;

  int ii = iran % 3;
  int jj = (iran/3) % 3;
  int kk = iran / 9;

  return lattice[i-1+ii][j-1+jj][k-1+kk];
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, update ghost values before computing propensity
------------------------------------------------------------------------- */

double AppPotts3d26n::site_propensity(int i, int j, int k, int full)
{
  if (full) site_update_ghosts(i,j,k);

  // loop over possible events
  // only consider spin flips to neighboring site values different from self

  int oldstate = lattice[i][j][k];
  int ii,jj,kk,newstate;
  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++) {
	newstate = lattice[ii][jj][kk];
	if (newstate == oldstate) continue;
	lattice[i][j][k] = newstate;
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
   if proc owns full domain, neighbor sites may be across PBC
   if only working on sector, ignore neighbor sites outside sector
------------------------------------------------------------------------- */

void AppPotts3d26n::site_event(int i, int j, int k, int full)
{
  int ii,jj,kk,isite,flag,sites[27];

  // pick one event from total propensity and set spin to that value

  double threshhold = random->uniform() * propensity[ijk2site[i][j][k]];

  int oldstate = lattice[i][j][k];
  int newstate;
  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++) {
	newstate = lattice[ii][jj][kk];
	if (newstate == oldstate) continue;
	lattice[i][j][k] = newstate;
	efinal = site_energy(i,j,k);
	if (efinal <= einitial) prob += 1.0;
	else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
	if (prob >= threshhold) goto done;
      }

 done:

  // compute propensity changes for self and neighbor sites

  int nsites = 0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++) {
	flag = 1;
	if (full) ijkpbc(ii,jj,kk);
	else if (ii < nx_sector_lo || ii > nx_sector_hi || 
		 jj < ny_sector_lo || jj > ny_sector_hi ||
		 kk < nz_sector_lo || kk > nz_sector_hi) flag = 0;
	if (flag) {
	  isite = ijk2site[ii][jj][kk];
	  sites[nsites++] = isite;
	  propensity[isite] = site_propensity(ii,jj,kk,full);
	}
      }
	
  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   update neighbors of site if neighbors are ghost cells
   called by site_propensity() when single proc owns entire domain
------------------------------------------------------------------------- */

void AppPotts3d26n::site_update_ghosts(int i, int j, int k)
{
  int ii,jj,kk;

  if (i == 1) {
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	lattice[i-1][jj][kk] = lattice[nx_local][jj][kk];
  }
  if (i == nx_local) {
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	lattice[i+1][jj][kk] = lattice[1][jj][kk];
  }
  if (j == 1) {
    for (ii = i-1; ii <= i+1; ii++)
      for (kk = k-1; kk <= k+1; kk++)
	lattice[ii][j-1][kk] = lattice[ii][ny_local][kk];
  }
  if (j == ny_local) {
    for (ii = i-1; ii <= i+1; ii++)
      for (kk = k-1; kk <= k+1; kk++)
	lattice[ii][j+1][kk] = lattice[ii][1][kk];
  }
  if (k == 1) {
    for (ii = i-1; ii <= i+1; ii++)
      for (jj = j-1; jj <= j+1; jj++)
	lattice[ii][jj][k-1] = lattice[ii][jj][nz_local];
  }
  if (k == nz_local) {
    for (ii = i-1; ii <= i+1; ii++)
      for (jj = j-1; jj <= j+1; jj++)
	lattice[ii][jj][k+1] = lattice[ii][jj][1];
  }
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppPotts3d26n::site_clear_mask(char ***mask, int i, int j, int k)
{
  int ii,jj,kk;
  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	mask[ii][jj][kk] = 0;
}
