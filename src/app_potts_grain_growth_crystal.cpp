/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.

   This app author:
   John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
------------------------------------------------------------------------- */

#include "app_potts_grain_growth_crystal.h"
#include "cubic_symmetries.h"
#include "disorientation.h"
#include "error.h"
#include "hcp_symmetries.h"
#include "math.h"
#include "quaternion.h"
#include "random_park.h"
#include "stdlib.h"
#include "string.h"
#include <cmath>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsGrainGrowthCrystal::AppPottsGrainGrowthCrystal(SPPARKS *spk, int narg,
                                                       char **arg)
    : AppPotts(spk, narg, arg), symmetries(), q0(nullptr), qx(nullptr),
      qy(nullptr), qz(nullptr), unique_neigh(nullptr) {
  // parse arguments for PottsNeighOnly class only, not children

  if (strcmp(style, "potts/graingrowth/crystal") != 0)
    return;

  if (narg != 3)
    error->all(FLERR, "Illegal 'potts/graingrowth/crystal command");

  if (strcmp(arg[2], "cubic") != 0) {
    symmetries = CUBIC::get_symmetries();
  } else if (strcmp(arg[2], "hcp") != 0) {
    symmetries = HCP::get_symmetries();
  } else {
    error->all(FLERR, "Illegal 'potts/graingrowth/crystal command; expected "
                      "'cubic' or 'hcp'");
  }

  nspins = atoi(arg[2]);
  if (nspins <= 0)
    error->all(FLERR, "Illegal app_style command");

  /* quaternion array for each site */
  ndouble = 4;
  // adds quaternion arrays
  recreate_arrays();
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsGrainGrowthCrystal::init_app() {
  delete[] sites;
  delete[] unique;
  delete[] unique_neigh;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];
  unique_neigh = new int[1 + maxneigh];

  dt_sweep = 1.0 / maxneigh;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    // Create random orientation as each site
    vector<double> uq = quaternion::generate_random_unit_quaternions(1);
    q0[i] = uq[0];
    qx[i] = uq[1];
    qy[i] = uq[2];
    qz[i] = uq[3];
    // check validity of spins
    if (spin[i] < 1 || spin[i] > nspins)
      flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall)
    error->all(FLERR, "One or more sites have invalid values");
}

void AppPottsGrainGrowthCrystal::grow_app() {
  spin = iarray[0];
  q0 = darray[0];
  qx = darray[1];
  qy = darray[2];
  qz = darray[3];
}

void AppPottsGrainGrowthCrystal::flip_site(int i, const SiteState &s) {
  spin[i] = s.spin;
  q0[i] = s.q[0];
  qx[i] = s.q[1];
  qy[i] = s.q[2];
  qz[i] = s.q[3];
}

double AppPottsGrainGrowthCrystal::site_energy(int i) {
  double energy = 0.0;
  double hi_angle = 15.0;
  vector<double> qi{q0[i], qx[i], qy[i], qz[i]};
  for (int j = 0; j < numneigh[i]; j++) {
    int nj = neighbor[i][j];
    vector<double> qj{q0[nj], qx[nj], qy[nj], qz[nj]};
    double di = disorientation::compute_disorientation(symmetries, qi, qj);
    double ratio = di / hi_angle;
    if (ratio <= 0)
      continue;
    else if (ratio < 1.0)
      energy += ratio * (1 - log(ratio));
    else
      energy += 1.0;
  }
  // Each site carries half the grain boundary energy
  // Neighbor sites carry the other half
  return 0.5 * energy;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsGrainGrowthCrystal::site_event_rejection(int i,
                                                      RandomPark *random) {
  int oldstate = spin[i];
  // Old state
  SiteState s0(oldstate, {q0[i], qx[i], qy[i], qz[i]});

  double einitial = site_energy(i);

  // events = spin flips to neighboring site different than self

  int j, m, value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    int nj = neighbor[i][j];
    value = spin[nj];
    if (value == spin[i])
      continue;
    for (m = 0; m < nevent; m++)
      // don't flip to a site which is same
      if (value == unique[m])
        break;
    // true if didn't match existing spin value in m loop;
    // therefore break and goto next neighbor j
    if (m < nevent)
      continue;
    // save spin value and associated neighbor
    unique[nevent] = value;
    unique_neigh[nevent] = nj;
    nevent += 1;
  }

  if (nevent == 0)
    return;
  int iran = (int)(nevent * random->uniform());
  if (iran >= nevent)
    iran = nevent - 1;
  // spin[i] = unique[iran];
  int neighran = unique_neigh[iran];
  SiteState s1(unique[iran],
               {q0[neighran], qx[neighran], qy[neighran], qz[neighran]});
  flip_site(i, s1);
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    // spin[i] = oldstate;
    flip_site(i, s0);
  } else if (random->uniform() > exp((einitial - efinal) * t_inverse)) {
    // spin[i] = oldstate;
    flip_site(i, s0);
  }

  if (spin[i] != oldstate)
    naccept++;

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5 * numneigh[i])
      mask[i] = 1;
    if (spin[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
        mask[neighbor[i][j]] = 0;
  }
}

double AppPottsGrainGrowthCrystal::site_propensity(int) {
  error->all(FLERR, "Illegal potts/graingrowth/crystal solver used.  KMC not "
                    "implemented.  Use rKMC.");
  return -1.0;
}
