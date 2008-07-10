/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_2d_24n.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppPotts2d24n::AppPotts2d24n(SPPARKS *spk, int narg, char **arg) : 
  AppPotts2d(spk,narg,arg)
{
  delevent = 0;
  delpropensity = 2;

  // parse any remaining arguments

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sample_argument") == 0) {
      iarg ++;
    } else {
      error->all("Illegal app_style potts/2d/24n command");
    }
  }

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nxlo,nxhi,nylo,nyhi,
			    "applattice2d:lattice");

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assignment is independent of # of procs

  if (init_style == RANDOM) {
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
  } else if (init_style == READ) read_spins(spinfile);
}

/* ---------------------------------------------------------------------- */

AppPotts2d24n::~AppPotts2d24n()
{
  memory->destroy_2d_T_array(lattice,nxlo,nylo);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts2d24n::site_energy(int i, int j)
{
  int isite = lattice[i][j];
  int eng = 0;

  if (isite != lattice[i-2][j-2]) eng++;
  if (isite != lattice[i-2][j-1]) eng++;
  if (isite != lattice[i-2][j]) eng++;
  if (isite != lattice[i-2][j+1]) eng++;
  if (isite != lattice[i-2][j+2]) eng++;

  if (isite != lattice[i-1][j-2]) eng++;
  if (isite != lattice[i-1][j-1]) eng++;
  if (isite != lattice[i-1][j]) eng++;
  if (isite != lattice[i-1][j+1]) eng++;
  if (isite != lattice[i-1][j+2]) eng++;

  if (isite != lattice[i][j-2]) eng++;
  if (isite != lattice[i][j-1]) eng++;
  if (isite != lattice[i][j+1]) eng++;
  if (isite != lattice[i][j+2]) eng++;

  if (isite != lattice[i+1][j-2]) eng++;
  if (isite != lattice[i+1][j-1]) eng++;
  if (isite != lattice[i+1][j]) eng++;
  if (isite != lattice[i+1][j+1]) eng++;
  if (isite != lattice[i+1][j+2]) eng++;

  if (isite != lattice[i+2][j-2]) eng++;
  if (isite != lattice[i+2][j-1]) eng++;
  if (isite != lattice[i+2][j]) eng++;
  if (isite != lattice[i+2][j+1]) eng++;
  if (isite != lattice[i+2][j+2]) eng++;

  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
   if site changes, unset mask of all neighbor sites with affected propensity
------------------------------------------------------------------------- */

void AppPotts2d24n::site_event_rejection(int i, int j, RandomPark *random)
{
  int oldstate = lattice[i][j];
  double einitial = site_energy(i,j);

  // event = random spin

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  lattice[i][j] = iran;
  double efinal = site_energy(i,j);

  // event = random neighbor spin

  //int iran = (int) (24.0*random->uniform());
  //if (iran == 0) return lattice[i-2][j-2];
  //else if (iran == 1) return lattice[i-2][j-1];
  //else if (iran == 2) return lattice[i-2][j];
  //else if (iran == 3) return lattice[i-2][j+1];
  //else if (iran == 4) return lattice[i-2][j+2];
  //else if (iran == 5) return lattice[i-1][j-2];
  //else if (iran == 6) return lattice[i-1][j-1];
  //else if (iran == 7) return lattice[i-1][j];
  //else if (iran == 8) return lattice[i-1][j+1];
  //else if (iran == 9) return lattice[i-1][j+2];
  //else if (iran == 10) return lattice[i][j-2];
  //else if (iran == 11) return lattice[i][j-1];
  //else if (iran == 12) return lattice[i][j+1];
  //else if (iran == 13) return lattice[i][j+2];
  //else if (iran == 14) return lattice[i+1][j-2];
  //else if (iran == 15) return lattice[i+1][j-1];
  //else if (iran == 16) return lattice[i+1][j];
  //else if (iran == 17) return lattice[i+1][j+1];
  //else if (iran == 18) return lattice[i+1][j+2];
  //else if (iran == 19) return lattice[i+2][j-2];
  //else if (iran == 20) return lattice[i+2][j-1];
  //else if (iran == 21) return lattice[i+2][j];
  //else if (iran == 22) return lattice[i+2][j+1];
  //else return lattice[i+2][j+2];

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
    if (einitial < 12.0) mask[i][j] = 1;
    if (lattice[i][j] != oldstate) {
      mask[i-2][j-2] = 0;
      mask[i-2][j-1] = 0;
      mask[i-2][j] = 0;
      mask[i-2][j+1] = 0;
      mask[i-2][j+2] = 0;

      mask[i-1][j-2] = 0;
      mask[i-1][j-1] = 0;
      mask[i-1][j] = 0;
      mask[i-1][j+1] = 0;
      mask[i-1][j+2] = 0;

      mask[i][j-2] = 0;
      mask[i][j-1] = 0;
      mask[i][j+1] = 0;
      mask[i][j+2] = 0;

      mask[i+1][j-2] = 0;
      mask[i+1][j-1] = 0;
      mask[i+1][j] = 0;
      mask[i+1][j+1] = 0;
      mask[i+1][j+2] = 0;

      mask[i+2][j-2] = 0;
      mask[i+2][j-1] = 0;
      mask[i+2][j] = 0;
      mask[i+2][j+1] = 0;
      mask[i+2][j+2] = 0;
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

double AppPotts2d24n::site_propensity(int i, int j)
{
  // possible events = spin flips to neighboring site different than self

  int sites[25];
  int oldstate = lattice[i][j];
  int nevent = 0;

  int ii,jj;
  for (ii = i-2; ii <= i+2; ii++)
    for (jj = j-2; jj <= j+2; jj++)
      add_unique(oldstate,nevent,sites,ii,jj);

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

void AppPotts2d24n::site_event(int i, int j, int full, RandomPark *random)
{
  // pick one event from total propensity

  double threshhold = random->uniform() * propensity[ij2site[i][j]];

  // possible events = spin flips to neighboring site different than self
  // find one event by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  int sites[25];
  int oldstate = lattice[i][j];
  int nevent = 0;

  int ii,jj;
  for (ii = i-2; ii <= i+2; ii++)
    for (jj = j-2; jj <= j+2; jj++)
      add_unique(oldstate,nevent,sites,ii,jj);

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

  int iloop,jloop,isite;

  int nsites = 0;

  for (iloop = i-2; iloop <= i+2; iloop++)
    for (jloop = j-2; jloop <= j+2; jloop++) {
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
