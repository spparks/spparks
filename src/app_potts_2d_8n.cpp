/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
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

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppPotts2d8n::AppPotts2d8n(SPK *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Invalid app_style potts/2d/8n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  seed = atoi(arg[4]);
  random = new RandomPark(seed);

  masklimit = 4.0;

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nx_local+2,ny_local+2,
			    "applattice2d:lattice");

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assigment is independent of # of procs

  int i,j,ii,jj,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      isite = random->irandom(nspins);
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local)
	lattice[ii][jj] = isite;
    }
  }

  // setup communicator for ghost sites

  comm = new CommLattice2d(spk);
  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth);
}

/* ---------------------------------------------------------------------- */

AppPotts2d8n::~AppPotts2d8n()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
  delete comm;
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
   randomly pick new state for site
------------------------------------------------------------------------- */

int AppPotts2d8n::site_pick_random(int i, int j, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

int AppPotts2d8n::site_pick_local(int i, int j, double ran)
{
  int iran = (int) (8*ran) + 1;
  if (iran > 8) iran = 8;

  if (iran == 1) return lattice[i-1][j-1];
  else if (iran == 2) return lattice[i-1][j];
  else if (iran == 3) return lattice[i-1][j+1];
  else if (iran == 4) return lattice[i][j-1];
  else if (iran == 5) return lattice[i][j+1];
  else if (iran == 6) return lattice[i+1][j-1];
  else if (iran == 7) return lattice[i+1][j];
  else return lattice[i+1][j+1];
}

/* ----------------------------------------------------------------------
   add this neighbor spin to set of possible new spins
------------------------------------------------------------------------- */

void AppPotts2d8n::survey_neighbor(const int& ik, const int& jk, int& ns, int spins[], int nspins[]) const {
  int *spnt = spins;
  bool Lfound;

  Lfound = false;
  while (spnt < spins+ns) {
    if (jk == *spnt++) {
      Lfound = true;
      break;
    }
  }

  if (Lfound) {
    // If found, increment counter 
    nspins[spnt-spins-1]++;
  } else {
    // If not found, create new survey entry
    spins[ns] = jk;
    nspins[ns] = 1;
    ns++;
  }

}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, update ghost values before computing propensity
------------------------------------------------------------------------- */

double AppPotts2d8n::site_propensity(int i, int j, int full)
{
  if (full) site_update_ghosts(i,j);

  // loop over possible events
  // only consider spin flips to neighboring site values different from self

  int oldstate = lattice[i][j];
  int ii,jj,newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  // Data for each possible new spin
  // nspins no longer used, as it is recalculated by site_energy(),
  // which is somewhat wasteful, but more general.
  int ns,spins[9],nspins[9];
  ns = 0;

  for (ii = i-1; ii <= i+1; ii++) {
    for (jj = j-1; jj <= j+1; jj++) {
      newstate = lattice[ii][jj];
      if (newstate == oldstate) continue;
      survey_neighbor(oldstate,newstate,ns,spins,nspins);
    }
  }

  // Use survey to compute overall propensity
  
  for (int is=0;is<ns;is++) {
    lattice[i][j] = spins[is];
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
   if proc owns full domain, neighbor sites may be across PBC
   if only working on sector, ignore neighbor sites outside sector
------------------------------------------------------------------------- */

void AppPotts2d8n::site_event(int i, int j, int full)
{
  int ii,jj,iloop,jloop,isite,flag,sites[9];

  // pick one event from total propensity and set spin to that value

  double threshhold = random->uniform() * propensity[ij2site[i][j]];

  int oldstate = lattice[i][j];
  int newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  // Data for each possible new spin
  // nspins no longer used, as it is recalculated by site_energy(),
  // which is somewhat wasteful, but more general.

  int ns, spins[8],nspins[8];
  ns = 0;

  for (ii = i-1; ii <= i+1; ii++) {
    for (jj = j-1; jj <= j+1; jj++) {
      newstate = lattice[ii][jj];
      if (newstate == oldstate) continue;
      survey_neighbor(oldstate,newstate,ns,spins,nspins);
    }
  }

  // Use survey to pick new spin
  
  for (int is=0;is<ns;is++) {
    lattice[i][j] = spins[is];
    efinal = site_energy(i,j);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) goto done;
  }

 done:

  // compute propensity changes for self and neighbor sites

  int nsites = 0;

  for (iloop = i-1; iloop <= i+1; iloop++)
    for (jloop = j-1; jloop <= j+1; jloop++) {
      ii = iloop; jj = jloop;
      flag = 1;
      if (full) ijpbc(ii,jj);
      else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	       jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
      if (flag) {
	isite = ij2site[ii][jj];
	sites[nsites++] = isite;
	propensity[isite] = site_propensity(ii,jj,full);
      }
    }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   update neighbors of site if neighbors are ghost cells
   called by site_propensity() when single proc owns entire domain
------------------------------------------------------------------------- */

void AppPotts2d8n::site_update_ghosts(int i, int j)
{
  if (i == 1) {
    lattice[i-1][j-1] = lattice[nx_local][j-1];
    lattice[i-1][j] = lattice[nx_local][j];
    lattice[i-1][j+1] = lattice[nx_local][j+1];
  }
  if (i == nx_local) {
    lattice[i+1][j-1] = lattice[1][j-1];
    lattice[i+1][j] = lattice[1][j];
    lattice[i+1][j+1] = lattice[1][j+1];
  }
  if (j == 1) {
    lattice[i-1][j-1] = lattice[i-1][ny_local];
    lattice[i][j-1] = lattice[i][ny_local];
    lattice[i+1][j-1] = lattice[i+1][ny_local];
  }
  if (j == ny_local) {
    lattice[i-1][j+1] = lattice[i-1][1];
    lattice[i][j+1] = lattice[i][1];
    lattice[i+1][j+1] = lattice[i+1][1];
  }
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppPotts2d8n::site_clear_mask(char **mask, int i, int j)
{
  mask[i-1][j-1] = 0;
  mask[i-1][j] = 0;
  mask[i-1][j+1] = 0;
  mask[i][j-1] = 0;
  mask[i][j] = 0;
  mask[i][j+1] = 0;
  mask[i+1][j-1] = 0;
  mask[i+1][j] = 0;
  mask[i+1][j+1] = 0;
}
