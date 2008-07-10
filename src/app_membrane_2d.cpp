/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_membrane_2d.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NONE,LIPID,FLUID,PROTEIN};

/* ---------------------------------------------------------------------- */

AppMembrane2d::AppMembrane2d(SPPARKS *spk, int narg, char **arg) : 
  AppLattice2d(spk,narg,arg)
{
  // parse arguments

  if (narg != 7) error->all("Illegal app_style membrane/2d command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  w01 = atof(arg[3]);
  w11 = atof(arg[4]);
  mu = atof(arg[5]);
  int seed = atoi(arg[6]);
  random = new RandomPark(seed);

  // define lattice and partition it across processors

  procs2lattice();
  memory->create_2d_T_array(lattice,nxlo,nxhi,nylo,nyhi,
			    "applattice2d:lattice");

  // initialize my portion of lattice to LIPID

  int i,j;
  for (i = 1; i <= nx_local; i++)
    for (j = 1; j <= ny_local; j++)
      lattice[i][j] = LIPID;

  // setup interaction energy matrix
  // w11 = fluid-fluid interaction
  // w01 = fluid-protein interaction

  interact[LIPID][LIPID] = interact[PROTEIN][PROTEIN] = 0.0;
  interact[LIPID][FLUID] = interact[FLUID][LIPID] = 0.0;
  interact[LIPID][PROTEIN] = interact[PROTEIN][LIPID] = 0.0;
  interact[FLUID][FLUID] = -w11;
  interact[FLUID][PROTEIN] = interact[PROTEIN][FLUID] = -w01;
}

/* ---------------------------------------------------------------------- */

AppMembrane2d::~AppMembrane2d()
{
  delete random;
  memory->destroy_2d_T_array(lattice);
}

/* ---------------------------------------------------------------------- */

void AppMembrane2d::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"inclusion") == 0) {
    if (narg != 3) error->all("Illegal inclusion command");
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
  } else error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppMembrane2d::site_energy(int i, int j)
{
  int isite = lattice[i][j];
  double eng = 0.0;
  eng += interact[isite][lattice[i-1][j]];
  eng += interact[isite][lattice[i+1][j]];
  eng += interact[isite][lattice[i][j-1]];
  eng += interact[isite][lattice[i][j+1]];
  return eng;
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
   if site cannot change, set mask
------------------------------------------------------------------------- */

void AppMembrane2d::site_event_rejection(int i, int j, RandomPark *random)
{
  int oldstate = lattice[i][j];
  double einitial = site_energy(i,j);

  // event = PROTEIN never changes, flip between LIPID and FLUID

  if (lattice[i][j] == PROTEIN) {
    if (Lmask) mask[i][j] = 1;
    return;
  }

  if (random->uniform() < 0.5) lattice[i][j] = LIPID;
  else lattice[i][j] = FLUID;
  double efinal = site_energy(i,j);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i][j] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i][j] = oldstate;
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
   propensity for one event based on einitial,efinal
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity = Boltzmann factor
------------------------------------------------------------------------- */

double AppMembrane2d::site_propensity(int i, int j)
{
  // only event is a LIPID/FLUID flip

  int oldstate = lattice[i][j];
  if (oldstate == PROTEIN) return 0.0;
  int newstate = LIPID;
  if (oldstate == LIPID) newstate = FLUID;

  // compute energy difference between initial and final state

  double einitial = site_energy(i,j);
  lattice[i][j] = newstate;
  double efinal = site_energy(i,j);
  lattice[i][j] = oldstate;

  if (oldstate == LIPID) efinal -= mu;
  else if (oldstate == FLUID) efinal += mu;

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

void AppMembrane2d::site_event(int i, int j, int full, RandomPark *random)
{
  // only event is a LIPID/FLUID flip

  if (lattice[i][j] == PROTEIN) return;
  if (lattice[i][j] == LIPID) lattice[i][j] = FLUID;
  else lattice[i][j] = LIPID;

  if (full) update_ghost_sites(i,j);

  // compute propensity changes for self and neighbor sites

  int ii,jj,isite,sites[5];

  int nsites = 0;

  ii = i; jj = j;
  isite = ij2site[ii][jj];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(ii,jj);

  ii = i-1; jj = j;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i+1; jj = j;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i; jj = j-1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  ii = i; jj = j+1;
  if (full) ijpbc(ii,jj);
  isite = ij2site[ii][jj];
  if (isite >= 0) {
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(ii,jj);
  }

  solve->update(nsites,sites,propensity);
}
