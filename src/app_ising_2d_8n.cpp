/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_ising_2d_8n.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppIsing2d8n::AppIsing2d8n(SPPARKS *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 4) error->all("Illegal app_style ising/2d/8n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  int seed = atoi(arg[3]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nxlo,nxhi,nylo,nyhi,
			    "applattice2d:lattice");

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

AppIsing2d8n::~AppIsing2d8n()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsing2d8n::site_energy(int i, int j)
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

void AppIsing2d8n::site_event_rejection(int i, int j, RandomPark *random)
{
  int oldstate = lattice[i][j];
  double einitial = site_energy(i,j);

  // event = random up or down spin

  if (random->uniform() < 0.5) lattice[i][j] = 1;
  else lattice[i][j] = 2;
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
   propensity for one event based on einitial,efinal
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity = Boltzmann factor
------------------------------------------------------------------------- */

double AppIsing2d8n::site_propensity(int i, int j)
{
  // event = spin flip

  int oldstate = lattice[i][j];
  int newstate = 1;
  if (oldstate == 1) newstate = 2;

  // compute energy difference between initial and final state

  double einitial = site_energy(i,j);
  lattice[i][j] = newstate;
  double efinal = site_energy(i,j);
  lattice[i][j] = oldstate;

  if (efinal <= einitial) return 1.0;
  else if (temperature == 0.0) return 0.0;
  else return exp((einitial-efinal)*t_inverse);
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, adjust neighbor sites indices for PBC
   if proc owns sector, ignore non-updated neighbors (isite < 0)
------------------------------------------------------------------------- */

void AppIsing2d8n::site_event(int i, int j, int full, RandomPark *random)
{
  // event = spin flip

  if (lattice[i][j] == 1) lattice[i][j] = 2;
  else lattice[i][j] = 1;

  if (full) update_ghost_sites(i,j);

  // compute propensity changes for self and neighbor sites

  int iloop,jloop,ii,jj,isite,sites[9];

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
