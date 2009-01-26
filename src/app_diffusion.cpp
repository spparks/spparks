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
#include "app_diffusion.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED};

/* ---------------------------------------------------------------------- */

AppDiffusion::AppDiffusion(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;
  allow_masking = 0;

  // parse arguments

  if (narg < 3) error->all("Illegal app_style command");

  double fraction = atof(arg[1]);
  int seed = atoi(arg[2]);
  random = new RandomPark(seed);

  options(narg-3,&arg[3]);

  // define lattice and partition it across processors
  // esites must be large enough for 2 sites and 1st/2nd nearest neighbors

  create_lattice();
  esites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  echeck = NULL;

  // initialize my portion of lattice
  // each site = VACANT or OCCUPIED with fraction OCCUPIED
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  if (infile) read_file();

  else {
    std::map<int,int> hash;
    for (int i = 0; i < nlocal; i++)
      hash.insert(std::pair<int,int> (id[i],i));
    std::map<int,int>::iterator loc;
    
    int isite;
    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      if (random->uniform() < fraction) isite = OCCUPIED;
      else isite = VACANT;
      loc = hash.find(iglobal);
      if (loc != hash.end()) lattice[loc->second] = isite;
    }
  }
}

/* ---------------------------------------------------------------------- */

AppDiffusion::~AppDiffusion()
{
  delete random;
  delete [] esites;
  delete [] echeck;
}

/* ---------------------------------------------------------------------- */

void AppDiffusion::init_app()
{
  delete [] echeck;
  echeck = new int[nlocal];
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusion::site_energy(int i)
{
  // energy only non-zero for OCCUPIED sites

  if (lattice[i] == VACANT) return 0.0;

  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
------------------------------------------------------------------------- */

void AppDiffusion::site_event_rejection(int i, RandomPark *random)
{
  // exchange with random neighbor if vacant
  // null bin extends to maxneigh

  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];
  if (lattice[j] == OCCUPIED) return;

  // accept or reject via Boltzmann criterion

  double einitial = site_energy(i) + site_energy(j);
  lattice[i] = VACANT;
  lattice[j] = OCCUPIED;
  double efinal = site_energy(i) + site_energy(j);

  if (efinal <= einitial) {
  } else if (temperature == 0) {
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }    
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusion::site_propensity(int i)
{
  int j;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // for each exchange
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  if (lattice[i] == VACANT) return 0.0;

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

  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusion::site_event(int i, class RandomPark *random)
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

  // compute propensity changes for self and swap site and their 1,2 neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

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

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
    esites[nsites++] = isite;
    echeck[isite] = 1;
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
