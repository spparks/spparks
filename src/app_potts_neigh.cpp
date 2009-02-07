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
#include "app_potts_neigh.h"
#include "random_park.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsNeigh::AppPottsNeigh(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  dt_sweep = 1.0/maxneigh;
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
   flip to random neighbor spin with null bin
   null bin extends to size maxneigh
------------------------------------------------------------------------- */

void AppPottsNeigh::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // events = spin flips to neighboring site different than self

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  int iran = (int) (maxneigh*random->uniform());
  if (iran >= nevent) return;
  lattice[i] = unique[iran];
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  if (lattice[i] != oldstate) naccept++;

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (lattice[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}
