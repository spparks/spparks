/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_ising_3d_6n.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppIsing3d6n::AppIsing3d6n(SPK *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
  // parse arguments

  if (narg != 5) error->all("Invalid app_style ising/3d/6n command");

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

AppIsing3d6n::~AppIsing3d6n()
{
  delete random;
  memory->destroy_3d_T_array(lattice);
  delete comm;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsing3d6n::site_energy(int i, int j, int k)
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

int AppIsing3d6n::site_pick_random(int i, int j, int k, double ran)
{
  int iran = (int) (2*ran) + 1;
  if (iran > 2) iran = 2;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
------------------------------------------------------------------------- */

int AppIsing3d6n::site_pick_local(int i, int j, int k, double ran)
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

double AppIsing3d6n::site_propensity(int i, int j, int k)
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

void AppIsing3d6n::site_event(int i, int j, int k)
{
  // only event is to flip the spin

  if (lattice[i][j][k] == 1) lattice[i][j][k] = 2;
  else lattice[i][j][k] = 1;

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

void AppIsing3d6n::site_update_ghost(int i, int j, int k)
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

void AppIsing3d6n::site_clear_mask(char ***mask, int i, int j, int k)
{
  mask[i][j][k] = 0;
  mask[i-1][j][k] = 0;
  mask[i+1][j][k] = 0;
  mask[i][j-1][k] = 0;
  mask[i][j+1][k] = 0;
  mask[i][j][k-1] = 0;
  mask[i][j][k+1] = 0;
}
