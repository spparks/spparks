/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_ising_3d_26n.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppIsing3d26n::AppIsing3d26n(SPPARKS *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Illegal app_style ising/3d/26n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nz_global = atoi(arg[3]);
  int seed = atoi(arg[4]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_3d_T_array(lattice,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
			    "app:lattice");

  // initialize my portion of lattice
  // each site = one of 2 spins
  // loop over global list so assigment is independent of # of procs

  int i,j,k,ii,jj,kk,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      for (k = 1; k <= nz_global; k++) {
	kk = k - nz_offset;
	isite = random->irandom(2);
	if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local &&
	    kk >= 1 && kk <= nz_local)
	  lattice[ii][jj][kk] = isite;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

AppIsing3d26n::~AppIsing3d26n()
{
  delete random;
  memory->destroy_3d_T_array(lattice);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsing3d26n::site_energy(int i, int j, int k)
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

void AppIsing3d26n::site_event_rejection(int i, int j, int k,
					RandomPark *random)
{
  int oldstate = lattice[i][j][k];
  double einitial = site_energy(i,j,k);

  // event = random up or down spin

  if (random->uniform() < 0.5) lattice[i][j][k] = 1;
  else lattice[i][j][k] = 2;
  double efinal = site_energy(i,j,k);

  // even = random neighbor spin

  //int iran = (int) (26*random->uniform()) + 1;
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
    if (einitial < 2.0) mask[i][j][k] = 1;
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
   propensity for one event based on einitial,efinal
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity = Boltzmann factor
------------------------------------------------------------------------- */

double AppIsing3d26n::site_propensity(int i, int j, int k)
{
  // only event is a spin flip

  int oldstate = lattice[i][j][k];
  int newstate = 1;
  if (oldstate == 1) newstate = 2;

  // compute energy difference between initial and final state

  double einitial = site_energy(i,j,k);
  lattice[i][j][k] = newstate;
  double efinal = site_energy(i,j,k);
  lattice[i][j][k] = oldstate;

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

void AppIsing3d26n::site_event(int i, int j, int k, int full,
			       RandomPark *random)
{
  // event = spin flip

  if (lattice[i][j][k] == 1) lattice[i][j][k] = 2;
  else lattice[i][j][k] = 1;

  if (full) update_ghost_sites(i,j,k);

  // compute propensity changes for self and neighbor sites

  int iloop,jloop,kloop,ii,jj,kk,isite,sites[27];

  int nsites = 0;

  for (iloop = i-1; iloop <= i+1; iloop++)
    for (jloop = j-1; jloop <= j+1; jloop++)
      for (kloop = k-1; kloop <= k+1; kloop++) {
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
