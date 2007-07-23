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
   pick new state for site randomly
------------------------------------------------------------------------- */

int AppPotts2d8n::site_pick_random(int i, int j, double ran)
{
  int iran = (int) (nspins*ran) + 1;
  if (iran > nspins) iran = nspins;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
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
   compute total propensity of site
   propensity based on einitial,efinal for each possible event
   no energy change = propensity of 1
   downhill energy change = propensity of 1
   uphill energy change = propensity via Boltzmann factor
------------------------------------------------------------------------- */

double AppPotts2d8n::site_propensity(int i, int j)
{
  site_update_ghost(i,j);

  // loop over possible events
  // only consider spin flips to neighboring site values different from self

  int oldstate = lattice[i][j];
  int ii,jj,newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++) {
      newstate = lattice[ii][jj];
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

void AppPotts2d8n::site_event(int i, int j)
{
  // pick one event from total propensity and set spin to that value

  double threshhold = random->uniform() * propensity[ij2site(i,j)];

  int oldstate = lattice[i][j];
  int ii,jj,newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++) {
      newstate = lattice[ii][jj];
      if (newstate == oldstate) continue;
      lattice[i][j] = newstate;
      efinal = site_energy(i,j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
      if (prob >= threshhold) goto done;
    }

 done:

  // reset propensity for self and neighbor sites

  int sites[9];
  int iii,jjj;
  int m = 0;

  for (iii = i-1; iii <= i+1; iii++)
    for (jjj = j-1; jjj <= j+1; jjj++) {
      ijpbc(iii,jjj,ii,jj);
      sites[m] = ij2site(ii,jj);
      propensity[sites[m++]] = site_propensity(ii,jj);
    }

  solve->update(9,sites,propensity);
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
------------------------------------------------------------------------- */

void AppPotts2d8n::site_event_sector(int i, int j)
{
  // pick one event from total propensity and set spin to that value
  // only update local sites; ghosts and pbcs are ignored

  double threshhold = random->uniform() * propensity[ij2site(i,j)];

  int oldstate = lattice[i][j];
  int ii,jj,newstate;
  double einitial = site_energy(i,j);
  double efinal;
  double prob = 0.0;

  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++) {
      newstate = lattice[ii][jj];
      if (newstate == oldstate) continue;
      lattice[i][j] = newstate;
      efinal = site_energy(i,j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
      if (prob >= threshhold) goto done;
    }

 done:

  // reset propensity for self and neighbor sites

  int sites[9];
  int m = 0;
  int isite;

  for (ii = MAX(1,i-1); ii <= MIN(i+1,nx_local); ii++)
    for (jj = MAX(1,j-1); jj <= MIN(j+1,ny_local); jj++) {
      isite = ij2site(ii,jj);
      sites[m] = isite;
      if (propensity[isite] >= 0.0) {
	propensity[isite] = site_propensity(ii,jj);
	m++;
      }
    }

  solve->update(m,sites,propensity);
}

/* ----------------------------------------------------------------------
  update neighbor cells of site that are global ghost cells
------------------------------------------------------------------------- */

void AppPotts2d8n::site_update_ghost(int i, int j)
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
