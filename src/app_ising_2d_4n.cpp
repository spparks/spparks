/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_ising_2d_4n.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

AppIsing2d4n::AppIsing2d4n(SPK *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 4) error->all("Invalid app_style ising/2d/4n command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  seed = atoi(arg[3]);
  random = new RandomPark(seed);

  masklimit = 2.0;

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nx_local+2,ny_local+2,
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

  // setup communicator for ghost sites

  comm = new CommLattice2d(spk);
  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth);
}

/* ---------------------------------------------------------------------- */

AppIsing2d4n::~AppIsing2d4n()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
  delete comm;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppIsing2d4n::site_energy(int i, int j)
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
   randomly pick new state for site
------------------------------------------------------------------------- */

int AppIsing2d4n::site_pick_random(int i, int j, double ran)
{
  int iran = (int) (2*ran) + 1;
  if (iran > 2) iran = 2;
  return iran;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

int AppIsing2d4n::site_pick_local(int i, int j, double ran)
{
  int iran = (int) (4*ran) + 1;
  if (iran > 4) iran = 4;

  if (iran == 1) return lattice[i-1][j];
  else if (iran == 2) return lattice[i+1][j];
  else if (iran == 3) return lattice[i][j-1];
  else return lattice[i][j+1];
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, update ghost values before computing propensity
------------------------------------------------------------------------- */

double AppIsing2d4n::site_propensity(int i, int j, int full)
{
  if (full) site_update_ghosts(i,j);

  // only event is a spin flip

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
   if proc owns full domain, neighbor sites may be across PBC
   if only working on sector, ignore neighbor sites outside sector
------------------------------------------------------------------------- */

void AppIsing2d4n::site_event(int i, int j, int full)
{
  int ii,jj,isite,flag,sites[5];

  // only event is a spin flip

  if (lattice[i][j] == 1) lattice[i][j] = 2;
  else lattice[i][j] = 1;

  // compute propensity changes for self and neighbor sites

  int nsites = 0;

  ii = i+1; jj = j;
  flag = 1;
  if (full) ijpbc(ii,jj);
  else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	   jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
  if (flag) {
    isite = ij2site[ii][jj];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,full);
  }

  ii = i+1; jj = j;
  flag = 1;
  if (full) ijpbc(ii,jj);
  else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	   jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
  if (flag) {
    isite = ij2site[ii][jj];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,full);
  }

  ii = i; jj = j;
  isite = ij2site[ii][jj];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(ii,jj,full);

  ii = 1; jj = j-1;
  flag = 1;
  if (full) ijpbc(ii,jj);
  else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	   jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
  if (flag) {
    isite = ij2site[ii][jj];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,full);
  }

  ii = 1; jj = j+1;
  flag = 1;
  if (full) ijpbc(ii,jj);
  else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	   jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
  if (flag) {
    isite = ij2site[ii][jj];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,full);
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
  update neighbors of site which has ghost cells for neighbors
  called bt site_propensity() when single proc owns entire domain
------------------------------------------------------------------------- */

void AppIsing2d4n::site_update_ghosts(int i, int j)
{
  if (i == 1) lattice[i-1][j] = lattice[nx_local][j];
  if (i == nx_local) lattice[i+1][j] = lattice[1][j];
  if (j == 1) lattice[i][j-1] = lattice[i][ny_local];
  if (j == ny_local) lattice[i][j+1] = lattice[i][1];
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppIsing2d4n::site_clear_mask(char **mask, int i, int j)
{
  mask[i][j] = 0;
  mask[i-1][j] = 0;
  mask[i+1][j] = 0;
  mask[i][j-1] = 0;
  mask[i][j+1] = 0;
}
