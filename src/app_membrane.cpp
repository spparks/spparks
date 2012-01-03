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
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NONE,LIPID,FLUID,PROTEIN};

/* ---------------------------------------------------------------------- */

AppMembrane::AppMembrane(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 0;
  delpropensity = 1;
  delevent = 0;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 0;
  numrandom = 1;
  dt_sweep = 1.0/2.0;

  create_arrays();

  // parse arguments

  if (narg != 4) error->all(FLERR,"Illegal app_style command");

  w01 = atof(arg[1]);
  w11 = atof(arg[2]);
  mu = atof(arg[3]);

  // setup interaction energy matrix
  // w11 = fluid-fluid interaction
  // w01 = fluid-protein interaction

  interact[LIPID][LIPID] = interact[PROTEIN][PROTEIN] = 0.0;
  interact[LIPID][FLUID] = interact[FLUID][LIPID] = 0.0;
  interact[LIPID][PROTEIN] = interact[PROTEIN][LIPID] = 0.0;
  interact[FLUID][FLUID] = -w11;
  interact[FLUID][PROTEIN] = interact[PROTEIN][FLUID] = -w01;

  sites = NULL;
}

/* ---------------------------------------------------------------------- */

AppMembrane::~AppMembrane()
{
  delete [] sites;
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppMembrane::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"inclusion") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal inclusion command");
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
      if (sqrt(rsq) < r) spin[i] = PROTEIN;
    }
  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppMembrane::grow_app()
{
  spin = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppMembrane::init_app()
{
  delete [] sites;
  sites = new int[1 + maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (spin[i] < LIPID || spin[i] > PROTEIN) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppMembrane::site_energy(int i)
{
  int isite = spin[i];
  double eng = 0.0;
  for (int j = 0; j < numneigh[i]; j++)
    eng += interact[isite][spin[neighbor[i][j]]];
  return eng;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   null bin extends to size 2
------------------------------------------------------------------------- */

void AppMembrane::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy(i);

  // event = PROTEIN never changes, flip between LIPID and FLUID

  if (spin[i] == PROTEIN) {
    if (Lmask) mask[i] = 1;
    return;
  }

  if (random->uniform() < 0.5) spin[i] = LIPID;
  else spin[i] = FLUID;
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    spin[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    spin[i] = oldstate;
  }

  if (spin[i] != oldstate) naccept++;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppMembrane::site_propensity(int i)
{
  // only event is a LIPID/FLUID flip

  int oldstate = spin[i];
  if (oldstate == PROTEIN) return 0.0;
  int newstate = LIPID;
  if (oldstate == LIPID) newstate = FLUID;

  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  double einitial = site_energy(i);
  spin[i] = newstate;
  double efinal = site_energy(i);
  spin[i] = oldstate;

  if (oldstate == LIPID) efinal -= mu;
  else if (oldstate == FLUID) efinal += mu;

  if (efinal <= einitial) return 1.0;
  else if (temperature == 0.0) return 0.0;
  else return exp((einitial-efinal)*t_inverse);
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppMembrane::site_event(int i, RandomPark *random)
{
  // only event is a LIPID/FLUID flip

  if (spin[i] == PROTEIN) return;
  if (spin[i] == LIPID) spin[i] = FLUID;
  else spin[i] = LIPID;

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
