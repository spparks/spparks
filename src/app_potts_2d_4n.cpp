/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_2d_4n.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppPotts2d4n::AppPotts2d4n(SPK *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Invalid app_style potts/2d/4n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nspins = atoi(arg[3]);
  seed = atoi(arg[4]);
  random = new RandomPark(seed);

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

AppPotts2d4n::~AppPotts2d4n()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
  delete comm;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts2d4n::site_energy(int i, int j)
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
   pick new state for site randomly
------------------------------------------------------------------------- */

int AppPotts2d4n::site_pick_random(int i, int j, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
------------------------------------------------------------------------- */

int AppPotts2d4n::site_pick_local(int i, int j, double ran)
{
  int iran = (int) (4*ran) + 1;
  if (iran > 4) iran = 4;

  if (iran == 1) return lattice[i-1][j];
  else if (iran == 2) return lattice[i+1][j];
  else if (iran == 3) return lattice[i][j-1];
  else return lattice[i][j+1];
}

/* ----------------------------------------------------------------------
   compute total propensity of site
   propensity based on einitial,efinal for each possible event
   no energy change = propensity of 1
   downhill energy change = propensity of 1
   uphill energy change = propensity via Boltzmann factor
------------------------------------------------------------------------- */

double AppPotts2d4n::site_propensity(int i, int j)
{
  site_update_ghost(i,j);

  // loop over possible events
  // only consider spin flips to neighboring site values different from self

  int oldstate = lattice[i][j];
  int newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < 4; m++) {
    if (m == 0) newstate = lattice[i-1][j];
    else if (m == 1) newstate = lattice[i+1][j];
    else if (m == 2) newstate = lattice[i][j-1];
    else newstate = lattice[i][j+1];
    if (newstate == oldstate) continue;
    lattice[i][j] = newstate;
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
------------------------------------------------------------------------- */

void AppPotts2d4n::site_event(int i, int j)
{
  // pick one event from total propensity and set spin to that value

  double threshhold = random->uniform() * propensity[ij2site(i,j)];

  int oldstate = lattice[i][j];
  int newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (int m = 0; m < 4; m++) {
    if (m == 0) newstate = lattice[i-1][j];
    else if (m == 1) newstate = lattice[i+1][j];
    else if (m == 2) newstate = lattice[i][j-1];
    else newstate = lattice[i][j+1];
    if (newstate == oldstate) continue;
    lattice[i][j] = newstate;
    efinal = site_energy(i,j);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  // reset propensity for self and neighbor sites

  int sites[5];
  int ii,jj;

  ijpbc(i-1,j,ii,jj);
  sites[0] = ij2site(ii,jj);
  propensity[sites[0]] = site_propensity(ii,jj);

  ijpbc(i+1,j,ii,jj);
  sites[1] = ij2site(ii,jj);
  propensity[sites[1]] = site_propensity(ii,jj);

  sites[2] = ij2site(i,j);
  propensity[sites[2]] = site_propensity(i,j);

  ijpbc(i,j-1,ii,jj);
  sites[3] = ij2site(ii,jj);
  propensity[sites[3]] = site_propensity(ii,jj);

  ijpbc(i,j+1,ii,jj);
  sites[4] = ij2site(ii,jj);
  propensity[sites[4]] = site_propensity(ii,jj);

  solve->update(5,sites,propensity);
}

/* ----------------------------------------------------------------------
  update neighbor cells of site that are global ghost cells
------------------------------------------------------------------------- */

void AppPotts2d4n::site_update_ghost(int i, int j)
{
  if (i == 1) lattice[i-1][j] = lattice[nx_local][j];
  if (i == nx_local) lattice[i+1][j] = lattice[1][j];
  if (j == 1) lattice[i][j-1] = lattice[i][ny_local];
  if (j == ny_local) lattice[i][j+1] = lattice[i][1];
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppPotts2d4n::site_clear_mask(char **mask, int i, int j)
{
  mask[i][j] = 0;
  mask[i-1][j] = 0;
  mask[i+1][j] = 0;
  mask[i][j-1] = 0;
  mask[i][j+1] = 0;
}
