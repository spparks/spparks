/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_membrane.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

enum{NONE,LIPID,FLUID,PROTEIN};

/* ---------------------------------------------------------------------- */

AppMembrane::AppMembrane(SPK *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 6) error->all("Invalid app_style membrane command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  w01 = atof(arg[3]);
  w11 = atof(arg[4]);
  seed = atoi(arg[5]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors
  
  procs2lattice();
  memory->create_2d_T_array(lattice,nx_local+2,ny_local+2,
			    "applattice2d:lattice");

  // initialize my portion of lattice
  // each site = LIPID
  // loop over global list so assigment is independent of # of procs

  int i,j,ii,jj,isite;
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    for (j = 1; j <= ny_global; j++) {
      jj = j - ny_offset;
      if (ii >= 1 && ii <= nx_local && jj >= 1 && jj <= ny_local)
	lattice[ii][jj] = LIPID;
    }
  }

  // setup communicator for ghost sites

  comm = new CommLattice2d(spk);
  comm->init(nx_local,ny_local,procwest,proceast,procsouth,procnorth);
}

/* ---------------------------------------------------------------------- */

AppMembrane::~AppMembrane()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
  delete comm;
}

/* ---------------------------------------------------------------------- */

void AppMembrane::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"inclusion") == 0) {
    if (narg != 3) error->all("Invalid inclusion command");
    int xc = atoi(arg[0]);
    int yc = atoi(arg[1]);
    double r = atof(arg[2]);
    int i,j,iglobal,jglobal;
    double rsq;
    for (i = 1; i <= nx_local; i++) {
      iglobal = i + nx_offset;
      for (j = 1; j <= ny_local; j++) {
	jglobal = j + ny_offset;
	rsq = (iglobal-xc)*(iglobal-xc) + (jglobal-yc)*(jglobal-yc);
	if (sqrt(rsq) < r) lattice[i][j] = PROTEIN;
      }
    }
  } else error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppMembrane::site_energy(int i, int j)
{
  int isite,jsite,itau,ieta,jtau,jeta;
  double eng = 0.0;

  isite = lattice[i][j];
  if (isite == FLUID) itau = 1;
  else itau = 0;
  if (isite == PROTEIN) ieta = 0;
  else ieta = 1;

  for (int m = 0; m < 8; m++) {
    //    if (m == 0) jsite = lattice[i-1][j];
    //    else if (m == 1) jsite = lattice[i+1][j];
    //    else if (m == 2) jsite = lattice[i][j-1];
    //    else jsite = lattice[i][j+1];

    if (m == 0) jsite = lattice[i-1][j-1];
    else if (m == 1) jsite = lattice[i-1][j];
    else if (m == 2) jsite = lattice[i-1][j+1];
    else if (m == 3) jsite = lattice[i][j-1];
    else if (m == 4) jsite = lattice[i][j+1];
    else if (m == 5) jsite = lattice[i+1][j-1];
    else if (m == 6) jsite = lattice[i+1][j];
    else jsite = lattice[i+1][j+1];

    if (jsite == FLUID) jtau = 1;
    else jtau = 0;
    if (jsite == PROTEIN) jeta = 0;
    else jeta = 1;

    // w11 = fluid-fluid interaction
    // w01 = fluid-protein interaction

    eng -= w11 * (itau*jtau*ieta*jeta);
    eng -= w01 * (itau*ieta*(1-jeta) + jtau*jeta*(1-ieta));
  }

  return eng;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site
------------------------------------------------------------------------- */

int AppMembrane::site_pick_random(int i, int j, double ran)
{
  if (lattice[i][j] == PROTEIN) return PROTEIN;
  if (ran < 0.5) return LIPID;
  else return FLUID;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

int AppMembrane::site_pick_local(int i, int j, double ran)
{
  if (lattice[i][j] == PROTEIN) return PROTEIN;
  if (ran < 0.5) return LIPID;
  else return FLUID;
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, update ghost values before computing propensity
------------------------------------------------------------------------- */

double AppMembrane::site_propensity(int i, int j, int full)
{
  if (full) site_update_ghosts(i,j);

  // only event is a LIPID/FLUID flip

  int oldstate = lattice[i][j];
  int newstate = LIPID;
  if (oldstate == LIPID) newstate = FLUID;

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

void AppMembrane::site_event(int i, int j, int full)
{
  int ii,jj,isite,flag,sites[5];

  // only event is a LIPID/FLUID flip

  if (lattice[i][j] == LIPID) lattice[i][j] = FLUID;
  else lattice[i][j] = LIPID;

  // compute propensity changes for self and neighbor sites

  int nsites = 0;

  ii = i; jj = j;
  isite = ij2site[ii][jj];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(ii,jj,full);

  ii = i-1; jj = j;
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

  ii = i; jj = j-1;
  flag = 1;
  if (full) ijpbc(ii,jj);
  else if (ii < nx_sector_lo || ii > nx_sector_hi || 
	   jj < ny_sector_lo || jj > ny_sector_hi) flag = 0;
  if (flag) {
    isite = ij2site[ii][jj];
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj,full);
  }

  ii = i; jj = j+1;
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
   update neighbors of site if neighbors are ghost cells
   called by site_propensity() when single proc owns entire domain
------------------------------------------------------------------------- */

void AppMembrane::site_update_ghosts(int i, int j)
{
  if (i == 1) lattice[i-1][j] = lattice[nx_local][j];
  if (i == nx_local) lattice[i+1][j] = lattice[1][j];
  if (j == 1) lattice[i][j-1] = lattice[i][ny_local];
  if (j == ny_local) lattice[i][j+1] = lattice[i][1];
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
------------------------------------------------------------------------- */

void AppMembrane::site_clear_mask(char **mask, int i, int j)
{
  mask[i][j] = 0;
  mask[i-1][j] = 0;
  mask[i+1][j] = 0;
  mask[i][j-1] = 0;
  mask[i][j+1] = 0;
}
