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
#include "string.h"
#include "stdlib.h"
#include "app_relax.h"
#include "potential.h"
#include "pair.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppRelax::AppRelax(SPPARKS *spk, int narg, char **arg) : 
  AppOffLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 0;
  allow_kmc = 1;
  allow_rejection = 1;
  dt_sweep = 1.0;

  create_arrays();

  // parse arguments

  if (narg != 2) error->all(FLERR,"Illegal app_style command");

  delta = atof(arg[1]);
  deltasq = delta*delta;
}

/* ---------------------------------------------------------------------- */

AppRelax::~AppRelax() {}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppRelax::grow_app()
{
  type = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppRelax::init_app()
{
  potential->init();
  pair = potential->pair;
  if (pair == NULL) error->all(FLERR,"App relax requires a pair potential");

  delpropensity = pair->cutoff;
  delevent = delta;

  int ntypes = potential->pair->ntypes;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (type[i] < 1 || type[i] > ntypes) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute energy of site
   form a neighbor list first
------------------------------------------------------------------------- */

double AppRelax::site_energy(int i)
{
  neighbor(i,pair->cutoff);
  return site_energy_neighbor(i);
}

/* ----------------------------------------------------------------------
   compute energy of site
   assume a neighbor list already exists
------------------------------------------------------------------------- */

double AppRelax::site_energy_neighbor(int i)
{
  return pair->energy(i,numneigh,neighs,xyz,type);
}

/* ----------------------------------------------------------------------
   perform a particle event with Metropolis algorithm
   event = translation move of up to delta
------------------------------------------------------------------------- */

void AppRelax::site_event_rejection(int i, RandomPark *random)
{
  double xold[3];
  double dx,dy,dz;

  double rc = pair->cutoff;
  neighbor(i,rc+delta);
  double einitial = site_energy_neighbor(i);

  xold[0] = xyz[i][0];
  xold[1] = xyz[i][1];
  xold[2] = xyz[i][2];

  double rsq = 1.0e20;
  while (rsq > deltasq) {
    dx = delta * (random->uniform() - 0.5);
    dy = delta * (random->uniform() - 0.5);
    if (dimension == 3) dz = delta * (random->uniform() - 0.5);
    else dz = 0.0;
    rsq = dx*dx + dy*dy + dz*dz;
  }
  xyz[i][0] += dx;
  xyz[i][1] += dy;
  xyz[i][2] += dz;

  double efinal = site_energy_neighbor(i);

  // accept or reject via Boltzmann criterion

  int success = 0;

  if (efinal <= einitial) {
    success = 1;
  } else if (temperature == 0.0) {
    xyz[i][0] = xold[0];
    xyz[i][1] = xold[1];
    xyz[i][2] = xold[2];
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    xyz[i][0] = xold[0];
    xyz[i][1] = xold[1];
    xyz[i][2] = xold[2];
  } else success = 1;

  if (success) {
    move(i);
    naccept++;
  }
}
