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
#include "app_potts_size.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsSize::AppPottsSize(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 1;
  allow_rejection = 0;
  allow_masking = 0;
  clustersizecol = 0;

  recreate_arrays();

  // parse arguments for PottsNeigh class only, not children

  if (strcmp(style,"potts/size") != 0) return;

  if (narg != 2) error->all(FLERR,"Illegal app_style command");
  
  nspins = atoi(arg[1]);
  if (nspins <= 0) error->all(FLERR,"Illegal app_style command");
  dt_sweep = 1.0/nspins;
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsSize::grow_app()
{
  spin = iarray[0];
  size = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsSize::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (spin[i] < 1 || spin[i] > nspins) flag = 1;
    if (size[i] < 0.0) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppPottsSize::site_propensity(int i)
{
  // events = spin flips to neighboring site different than self
  // disallow wild flips = flips to value different than all neighs

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == spin[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  // for each flip:
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor
  // TODO: cluster size scaling not defined yet

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double efinal,scale;
  double prob = 0.0;

  for (m = 0; m < nevent; m++) {
    spin[i] = unique[m];
    efinal = site_energy(i);
    if (efinal <= einitial) {
      //      scale = 1.0 + size[i];
      //      prob += 1.0/scale;
      prob += 1.0;
    } else if (temperature > 0.0) {
      // scale = 1.0 + size[i];
      // prob += exp((einitial-efinal)*(t_inverse/scale));
      prob += exp((einitial-efinal)*(t_inverse));
    }
  }

  spin[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPottsSize::site_event(int i, RandomPark *random)
{
  int j,m,value;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double efinal,scale;

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double prob = 0.0;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == oldstate) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
    
  // TODO: cluster size scaling not defined yet

    spin[i] = value;
    efinal = site_energy(i);
    if (efinal <= einitial) {
      // scale = 1.0 + size[i];
      // prob += 1.0/scale;
      prob += 1.0;
    } else if (temperature > 0.0) {
      // scale = 1.0 + size[i];
      // prob += exp((einitial-efinal)*(t_inverse/scale));
      prob += exp((einitial-efinal)*(t_inverse));
    }
    if (prob >= threshhold) break;
  }

  // compute propensity changes for self and neighbor sites
  // ignore update of neighbor sites with isite < 0

  int nsites = 0;
  int isite = i2site[i];
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

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
------------------------------------------------------------------------- */

void *AppPottsSize::extract_app(char *name)
{
  if (strcmp(name,"nspins") == 0) return (void *) &nspins;
  if (strcmp(name,"size") == 0) return (void *) size ;
  return NULL;
}
