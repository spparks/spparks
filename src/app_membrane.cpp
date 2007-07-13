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
   pick new state for site randomly
------------------------------------------------------------------------- */

int AppMembrane::site_pick_random(int i, int j, double ran)
{
  if (lattice[i][j] == PROTEIN) return PROTEIN;
  if (ran < 0.5) return LIPID;
  else return FLUID;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values
------------------------------------------------------------------------- */

int AppMembrane::site_pick_local(int i, int j, double ran)
{
  if (lattice[i][j] == PROTEIN) return PROTEIN;
  if (ran < 0.5) return LIPID;
  else return FLUID;
}

/* ----------------------------------------------------------------------
   compute total propensity of site
   propensity based on einitial,efinal for each possible event
   no energy change = propensity of 1
   downhill energy change = propensity of 1
   uphill energy change = propensity via Boltzmann factor
------------------------------------------------------------------------- */

double AppMembrane::site_propensity(int i, int j)
{
  site_update_ghost(i,j);

  // only event is to flip the spin

  int oldstate = lattice[i][j];
  int newstate = 1;
  if (oldstate == 1) newstate = 2;

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
------------------------------------------------------------------------- */

void AppMembrane::site_event(int i, int j)
{
  // only event is to flip the spin

  if (lattice[i][j] == 1) lattice[i][j] = 2;
  else lattice[i][j] = 1;

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

void AppMembrane::site_update_ghost(int i, int j)
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
