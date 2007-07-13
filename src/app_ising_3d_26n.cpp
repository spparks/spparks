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

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppIsing3d26n::AppIsing3d26n(SPK *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Invalid app_style ising/3d/26n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nz_global = atoi(arg[3]);
  seed = atoi(arg[4]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_3d_T_array(lattice,nx_local+2,ny_local+2,nz_local+2,
			    "applattice3d:lattice");

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

  // setup communicator for ghost sites

  comm = new CommLattice3d(spk);
  comm->init(nx_local,ny_local,nz_local,
	     procwest,proceast,procsouth,procnorth,procdown,procup);
}

/* ---------------------------------------------------------------------- */

AppIsing3d26n::~AppIsing3d26n()
{
  delete random;
  memory->destroy_3d_T_array(lattice);
  delete comm;
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
   pick new state for site randomly
------------------------------------------------------------------------- */

int AppIsing3d26n::site_pick_random(int i, int j, int k, double ran)
{
  int iran = (int) (2*ran) + 1;
  if (iran > 2) iran = 2;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
------------------------------------------------------------------------- */

int AppIsing3d26n::site_pick_local(int i, int j, int k, double ran)
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
   compute total propensity of site
   propensity based on einitial,efinal for each possible event
   no energy change = propensity of 1
   downhill energy change = propensity of 1
   uphill energy change = propensity via Boltzmann factor
------------------------------------------------------------------------- */

double AppIsing3d26n::site_propensity(int i, int j, int k)
{
  site_update_ghost(i,j,k);

  // only event is to flip the spin

  int oldstate = lattice[i][j][k];
  int newstate = 1;
  if (oldstate == 1) newstate = 2;

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
------------------------------------------------------------------------- */

void AppIsing3d26n::site_event(int i, int j, int k)
{
  // only event is to flip the sping

  if (lattice[i][j][k] == 1) lattice[i][j][k] = 2;
  else lattice[i][j][k] = 1;

  // reset propensity for self and neighbor sites

  int sites[27];
  int ii,jj,kk,iii,jjj,kkk;
  int m = 0;

  for (iii = i-1; iii <= i+1; iii++)
    for (jjj = j-1; jjj <= j+1; jjj++)
      for (kkk = k-1; kkk <= k+1; kkk++) {
	ijkpbc(iii,jjj,kkk,ii,jj,kk);
	sites[m] = ijk2site(ii,jj,kk);
	propensity[sites[m++]] = site_propensity(ii,jj,kk);
      }

  solve->update(27,sites,propensity);
}

/* ----------------------------------------------------------------------
  update neighbor cells of site that are global ghost cells
------------------------------------------------------------------------- */

void AppIsing3d26n::site_update_ghost(int i, int j, int k)
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

void AppIsing3d26n::site_clear_mask(char ***mask, int i, int j, int k)
{
  int ii,jj,kk;
  for (ii = i-1; ii <= i+1; ii++)
    for (jj = j-1; jj <= j+1; jj++)
      for (kk = k-1; kk <= k+1; kk++)
	mask[ii][jj][kk] = 0;
}
