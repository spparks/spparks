/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_membrane.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NONE,LIPID,FLUID,PROTEIN};

/* ---------------------------------------------------------------------- */

AppMembrane::AppMembrane(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // parse arguments

  if (narg < 5) error->all("Illegal app_style membrane command");

  w01 = atof(arg[1]);
  w11 = atof(arg[2]);
  mu = atof(arg[3]);
  seed = atoi(arg[4]);
  random = new RandomPark(seed);

  options(narg-5,&arg[5]);

  // setup interaction energy matrix
  // w11 = fluid-fluid interaction
  // w01 = fluid-protein interaction

  interact[LIPID][LIPID] = interact[PROTEIN][PROTEIN] = 0.0;
  interact[LIPID][FLUID] = interact[FLUID][LIPID] = 0.0;
  interact[LIPID][PROTEIN] = interact[PROTEIN][LIPID] = 0.0;
  interact[FLUID][FLUID] = -w11;
  interact[FLUID][PROTEIN] = interact[PROTEIN][FLUID] = -w01;

  // define lattice and partition it across processors

  create_lattice();
  sites = new int[1 + maxneigh];

  // initialize my sites to LIPID

  for (int i = 0; i < nlocal; i++) lattice[i] = LIPID;
}

/* ---------------------------------------------------------------------- */

AppMembrane::~AppMembrane()
{
  delete random;
  delete [] sites;
  delete comm;
}

/* ---------------------------------------------------------------------- */

void AppMembrane::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"inclusion") == 0) {
    if (narg != 4) error->all("Illegal inclusion command");
    double xc = atof(arg[0]);
    double yc = atof(arg[1]);
    double zc = atof(arg[2]);
    double r = atof(arg[3]);

    double dx,dy,dz,rsq;
    for (int i = 0; i < nlocal; i++) {
      dx = xyz[i][0] - xc;
      dy = xyz[i][1] - yc;
      dz = xyz[i][2] - zc;
      rsq = dx*dx + dy*dy + dz*dz;
      if (sqrt(rsq) < r) lattice[i] = PROTEIN;
    }
  } else error->all("Unrecognized command");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppMembrane::site_energy(int i)
{
  int isite = lattice[i];
  double eng = 0.0;
  for (int j = 0; j < numneigh[i]; j++)
    eng += interact[isite][lattice[neighbor[i][j]]];
  return eng;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site
------------------------------------------------------------------------- */

void AppMembrane::site_pick_random(int i, double ran)
{
  if (lattice[i] == PROTEIN) return;
  if (ran < 0.5) lattice[i] = LIPID;
  else lattice[i] = FLUID;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

void AppMembrane::site_pick_local(int i, double ran)
{
  if (lattice[i] == PROTEIN) return;
  if (ran < 0.5) lattice[i] = LIPID;
  else lattice[i] = FLUID;
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, there are no ghosts, so ignore full flag
------------------------------------------------------------------------- */

double AppMembrane::site_propensity(int i, int full)
{
  // only event is a LIPID/FLUID flip

  int oldstate = lattice[i];
  if (oldstate == PROTEIN) return 0.0;
  int newstate = LIPID;
  if (oldstate == LIPID) newstate = FLUID;

  double einitial = site_energy(i);
  lattice[i] = newstate;
  double efinal = site_energy(i);
  lattice[i] = oldstate;

  if (oldstate == LIPID) efinal -= mu;
  else if (oldstate == FLUID) efinal += mu;
  double delta = einitial - efinal;

  if (efinal <= einitial) return 1.0;
  else if (temperature == 0.0) return 0.0;
  else return exp(delta*t_inverse);
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, there are no ghosts, so ignore full flag
   if proc owns sector, ignore neighbor sites that are ghosts
------------------------------------------------------------------------- */

void AppMembrane::site_event(int i, int full)
{
  int j,m,isite;

  // only event is a LIPID/FLUID flip

  if (lattice[i] == PROTEIN) return;
  if (lattice[i] == LIPID) lattice[i] = FLUID;
  else lattice[i] = LIPID;

  // compute propensity changes for self and neighbor sites

  int nsites = 0;
  isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i,full);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m,full);
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
  OK to clear ghost site mask values
------------------------------------------------------------------------- */

void AppMembrane::site_clear_mask(char *mask, int i)
{
  mask[i] = 0;
  for (int j = 0; j < numneigh[i]; j++) mask[neighbor[i][j]] = 0;
}
