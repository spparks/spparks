/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_3d_6n.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppPotts3d6n::AppPotts3d6n(SPK *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
  // parse arguments

  if (narg != 6) error->all("Invalid app_style potts/3d/6n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nz_global = atoi(arg[3]);
  nspins = atoi(arg[4]);
  seed = atoi(arg[5]);
  random = new RandomPark(seed);

  masklimit = 3.0;

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

AppPotts3d6n::~AppPotts3d6n()
{
  delete random;
  memory->destroy_3d_T_array(lattice);
  delete comm;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts3d6n::site_energy(int i, int j, int k)
{
  int isite = lattice[i][j][k];
  int eng = 0;
  if (isite != lattice[i][j][k-1]) eng++;
  if (isite != lattice[i][j][k+1]) eng++;
  if (isite != lattice[i][j-1][k]) eng++;
  if (isite != lattice[i][j+1][k]) eng++;
  if (isite != lattice[i-1][j][k]) eng++;
  if (isite != lattice[i+1][j][k]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly
------------------------------------------------------------------------- */

int AppPotts3d6n::site_pick_random(int i, int j, int k, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
------------------------------------------------------------------------- */

int AppPotts3d6n::site_pick_local(int i, int j, int k, double ran)
{
  int iran = (int) (6*ran) + 1;
  if (iran > 6) iran = 6;

  if (iran == 1) return lattice[i-1][j][k];
  else if (iran == 2) return lattice[i+1][j][k];
  else if (iran == 3) return lattice[i][j-1][k];
  else if (iran == 4) return lattice[i][j+1][k];
  else if (iran == 5) return lattice[i][j][k-1];
  else return lattice[i][j][k+1];
}

/* ----------------------------------------------------------------------
   compute total propensity of site
   propensity based on einitial,efinal for each possible event
   no energy change = propensity of 1
   downhill energy change = propensity of 1
   uphill energy change = propensity via Boltzmann factor
------------------------------------------------------------------------- */

double AppPotts3d6n::site_propensity(int i, int j, int k)
{
  site_update_ghost(i,j,k);

  // loop over possible events
  // only consider spin flips to neighboring site values different from self

  int oldstate = lattice[i][j][k];
  int newstate;
  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < 6; m++) {
    if (m == 0) newstate = lattice[i-1][j][k];
    else if (m == 1) newstate = lattice[i+1][j][k];
    else if (m == 2) newstate = lattice[i][j-1][k];
    else if (m == 3) newstate = lattice[i][j+1][k];
    else if (m == 4) newstate = lattice[i][j][k-1];
    else newstate = lattice[i][j][k+1];
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
------------------------------------------------------------------------- */

void AppPotts3d6n::site_event(int i, int j, int k)
{
  // pick one event from total propensity and set spin to that value

  double threshhold = random->uniform() * propensity[ijk2site(i,j,k)];

  int oldstate = lattice[i][j][k];
  int newstate;
  double einitial = site_energy(i,j,k);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < 6; m++) {
    if (m == 0) newstate = lattice[i-1][j][k];
    else if (m == 1) newstate = lattice[i+1][j][k];
    else if (m == 2) newstate = lattice[i][j-1][k];
    else if (m == 3) newstate = lattice[i][j+1][k];
    else if (m == 4) newstate = lattice[i][j][k-1];
    else newstate = lattice[i][j][k+1];
    if (newstate == oldstate) continue;
    lattice[i][j][k] = newstate;
    efinal = site_energy(i,j,k);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  // reset propensity for self and neighbor sites

  int sites[7];
  int ii,jj,kk;

  ijkpbc(i-1,j,k,ii,jj,kk);
  sites[0] = ijk2site(ii,jj,kk);
  propensity[sites[0]] = site_propensity(ii,jj,kk);

  ijkpbc(i+1,j,k,ii,jj,kk);
  sites[1] = ijk2site(ii,jj,kk);
  propensity[sites[1]] = site_propensity(ii,jj,kk);

  ijkpbc(i,j-1,k,ii,jj,kk);
  sites[2] = ijk2site(ii,jj,kk);
  propensity[sites[2]] = site_propensity(ii,jj,kk);

  ijkpbc(i,j+1,k,ii,jj,kk);
  sites[3] = ijk2site(ii,jj,kk);
  propensity[sites[3]] = site_propensity(ii,jj,kk);

  ijkpbc(i,j,k-1,ii,jj,kk);
  sites[4] = ijk2site(ii,jj,kk);
  propensity[sites[4]] = site_propensity(ii,jj,kk);

  ijkpbc(i,j,k+1,ii,jj,kk);
  sites[5] = ijk2site(ii,jj,kk);
  propensity[sites[5]] = site_propensity(ii,jj,kk);

  sites[6] = ijk2site(i,j,k);
  propensity[sites[6]] = site_propensity(i,j,k);

  solve->update(7,sites,propensity);
}

/* ----------------------------------------------------------------------
  update neighbor cells of site that are global ghost cells
------------------------------------------------------------------------- */

void AppPotts3d6n::site_update_ghost(int i, int j, int k)
{
  if (i == 1) lattice[i-1][j][k] = lattice[nx_local][j][k];
  if (i == nx_local) lattice[i+1][j][k] = lattice[1][j][k];
  if (j == 1) lattice[i][j-1][k] = lattice[i][ny_local][k];
  if (j == ny_local) lattice[i][j+1][k] = lattice[i][1][k];
  if (k == 1) lattice[i][j][k-1] = lattice[i][j][nz_local];
  if (k == nz_local) lattice[i][j][k+1] = lattice[i][j][1];
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppPotts3d6n::site_clear_mask(char ***mask, int i, int j, int k)
{
  mask[i][j][k] = 0;
  mask[i-1][j][k] = 0;
  mask[i+1][j][k] = 0;
  mask[i][j-1][k] = 0;
  mask[i][j+1][k] = 0;
  mask[i][j][k-1] = 0;
  mask[i][j][k+1] = 0;
}
