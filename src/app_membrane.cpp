/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
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
  dt_sweep = 1.0/2.0;

  // parse arguments

  if (narg < 4) error->all("Illegal app_style command");

  w01 = atof(arg[1]);
  w11 = atof(arg[2]);
  mu = atof(arg[3]);

  options(narg-4,&arg[4]);

  // define lattice and partition it across processors

  create_lattice();
  sites = new int[1 + maxneigh];

  // setup interaction energy matrix
  // w11 = fluid-fluid interaction
  // w01 = fluid-protein interaction

  interact[LIPID][LIPID] = interact[PROTEIN][PROTEIN] = 0.0;
  interact[LIPID][FLUID] = interact[FLUID][LIPID] = 0.0;
  interact[LIPID][PROTEIN] = interact[PROTEIN][LIPID] = 0.0;
  interact[FLUID][FLUID] = -w11;
  interact[FLUID][PROTEIN] = interact[PROTEIN][FLUID] = -w01;

  // initialize my portion of lattice to LIPID

  if (infile) read_file();

  else {
    for (int i = 0; i < nlocal; i++) lattice[i] = LIPID;
  }
}

/* ---------------------------------------------------------------------- */

AppMembrane::~AppMembrane()
{
  delete [] sites;
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
   perform a site event with null bin rejection
   null bin extends to size 2
------------------------------------------------------------------------- */

void AppMembrane::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // event = PROTEIN never changes, flip between LIPID and FLUID

  if (lattice[i] == PROTEIN) {
    if (Lmask) mask[i] = 1;
    return;
  }

  if (random->uniform() < 0.5) lattice[i] = LIPID;
  else lattice[i] = FLUID;
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  if (lattice[i] != oldstate) naccept++;
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppMembrane::site_propensity(int i)
{
  // only event is a LIPID/FLUID flip

  int oldstate = lattice[i];
  if (oldstate == PROTEIN) return 0.0;
  int newstate = LIPID;
  if (oldstate == LIPID) newstate = FLUID;

  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  double einitial = site_energy(i);
  lattice[i] = newstate;
  double efinal = site_energy(i);
  lattice[i] = oldstate;

  if (oldstate == LIPID) efinal -= mu;
  else if (oldstate == FLUID) efinal += mu;

  if (efinal <= einitial) return 1.0;
  else if (temperature == 0.0) return 0.0;
  else return exp((einitial-efinal)*t_inverse);
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppMembrane::site_event(int i, RandomPark *random)
{
  // only event is a LIPID/FLUID flip

  if (lattice[i] == PROTEIN) return;
  if (lattice[i] == LIPID) lattice[i] = FLUID;
  else lattice[i] = LIPID;

  // compute propensity changes for self and neighbor sites
  // ignore update of neighbor sites with isite < 0

  int j,m,isite;

  int nsites = 0;
  isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
  }

  solve->update(nsites,sites,propensity);
}
