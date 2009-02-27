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
#include "app_diffusion_deposit.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED};

/* ---------------------------------------------------------------------- */

AppDiffusionDeposit::AppDiffusionDeposit(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;
  allow_rejection = 0;
  allow_masking = 0;

  // parse arguments

  if (narg < 3) error->all("Illegal app_style command");

  double surflevel = atof(arg[1]);
  deprate = atof(arg[2]);

  hopthresh = maxneigh/2.0 + 0.1;

  options(narg-3,&arg[3]);

  // don't allow parallel for now

  if (nprocs > 1)
    error->all("Cannot run app_style diffusion/deposit in parallel for now");

  // define lattice and partition it across processors
  // esites must be large enough for 2 sites and 1st/2nd nearest neighbors

  create_lattice();
  esites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  echeck = NULL;

  // initialize my portion of lattice
  // each site = VACANT or OCCUPIED with nlayers in x or xy OCCUPIED
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  RandomPark *random = new RandomPark(ranmaster->uniform());

  if (infile) read_file();

  else {
    std::map<int,int> hash;
    for (int i = 0; i < nlocal; i++)
      hash.insert(std::pair<int,int> (id[i],i));
    std::map<int,int>::iterator loc;
    
    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      loc = hash.find(iglobal);
      if (loc != hash.end()) {
	if (dimension == 2) {
	  if (xyz[loc->second][1] <= surflevel)
	    lattice[loc->second] = OCCUPIED;
	  else lattice[loc->second] = VACANT;
	} else {
	  if (xyz[loc->second][2] <= surflevel)
	    lattice[loc->second] = OCCUPIED;
	  else lattice[loc->second] = VACANT;
	}
      }
    }
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

AppDiffusionDeposit::~AppDiffusionDeposit()
{
  delete [] esites;
  delete [] echeck;
}

/* ---------------------------------------------------------------------- */

void AppDiffusionDeposit::init_app()
{
  delete [] echeck;
  echeck = new int[nlocal];
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusionDeposit::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != lattice[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusionDeposit::site_propensity(int i)
{
  int j;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // for each exchange
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  if (lattice[i] == VACANT) return 0.0;
  if (site_energy(i) > hopthresh) return 0.0;

  double einitial,efinal;
  double prob = 0.0;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANT) {
      einitial = site_energy(i) + site_energy(j);

      lattice[i] = VACANT;
      lattice[j] = OCCUPIED;
      efinal = site_energy(i) + site_energy(j);

      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);

      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // add in single deposition event, stored by site 0

  if (i == 0) prob += deprate;

  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusionDeposit::site_event(int i, class RandomPark *random)
{
  int j,jj,k,kk,m,mm;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double prob = 0.0;
  double einitial,efinal;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANT) {
      einitial = site_energy(i) + site_energy(j);
      lattice[i] = VACANT;
      lattice[j] = OCCUPIED;
      efinal = site_energy(i) + site_energy(j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
      if (prob >= threshhold) break;
      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // deposition event
  // find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (i == 0 && prob < threshhold) {
    if (dimension == 2) m = find_deposit_site_2d(random);
    else m = find_deposit_site_3d(random);
    lattice[m] = OCCUPIED;
    i = j = m;
  }

  // compute propensity changes for self and swap site and their 1,2 neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  for (k = 0; k < numneigh[j]; k++) {
    m = neighbor[j][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	esites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom for 3d lattices
------------------------------------------------------------------------- */

int AppDiffusionDeposit::find_deposit_site_3d(RandomPark *random)
{
  // pick a random position at top of box

  double xstart = boxxlo + (boxxhi-boxxlo)*random->uniform();
  double ystart = boxylo + (boxyhi-boxylo)*random->uniform();
  double zstart = boxzhi;

  // pick a random downward direction, make it a unit vector

  double dir[3];
  dir[0] = 2.0*random->uniform() - 1.0;
  dir[1] = 2.0*random->uniform() - 1.0;
  dir[2] = -0.5 + 0.5*random->uniform();

  double len = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  dir[0] /= len;
  dir[1] /= len;
  dir[2] /= len;

  // for each vacant site:
  // check if perp distance to incident line < lattice spacing
  //   account for periodic boundaries in XY
  // discard site if has no neighbors
  // neighbors across z boundary don't count
  // find site closest to start point that meets these criteria

  int closesite = -1;
  double closedistsq = 1.0e20;
  double dot,distsq,dx,dy,dz;
  int ncount,m,n;
  double delta[3],projection[3],offset[3];

  for (m = 0; m < nlocal; m++) {
    if (lattice[m] == OCCUPIED) continue;

    delta[0] = xyz[m][0] - xstart;
    delta[1] = xyz[m][1] - ystart;
    delta[2] = xyz[m][2] - zstart;

    dot = dir[0]*delta[0] + dir[1]*delta[1] + dir[2]*delta[2];
    projection[0] = xstart + dot*dir[0];
    projection[1] = ystart + dot*dir[1];
    projection[2] = zstart + dot*dir[2];

    offset[0] = xyz[m][0] - projection[0];
    offset[1] = xyz[m][1] - projection[1];
    offset[2] = xyz[m][2] - projection[2];

    distsq = offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
    if (distsq > nnspacingsq) continue;

    ncount = 0;
    for (n = 0; n < numneigh[m]; n++)
      if (neighbor[m][n] == OCCUPIED && 
	  xyz[neighbor[m][n]][2] != 0.0) ncount++;
    if (ncount == 0) continue;

    distsq = dx*dx + dy*dy + dz*dz;

    if (distsq < closedistsq) {
      closedistsq = distsq;
      closesite = m;
    }
  }

  return closesite;
}

/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom for 2d lattices
------------------------------------------------------------------------- */

int AppDiffusionDeposit::find_deposit_site_2d(RandomPark *random)
{
  return 0;
}
