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

/* ----------------------------------------------------------------------
   Contributing author: John Mitchell (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_strain_pin.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include <map>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsStrainPin::AppPottsStrainPin(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 1;
  allow_rejection = 0;
  allow_masking = 0;

  recreate_arrays();

  // parse arguments for PottsNeigh class only, not children

  if (strcmp(style,"potts/strain/with/pinning") != 0) return;

  if (narg != 2) error->all(FLERR,"Illegal app_style command");
  
  nspins = atoi(arg[1]);
  if (nspins <= 0) error->all(FLERR,"Illegal app_style command");
  dt_sweep = 1.0/nspins;
}

/* ----------------------------------------------------------------------
   input script for 'pin' command
------------------------------------------------------------------------- */

void AppPottsStrainPin::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"pin") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal pin command");
    pfraction = atof(arg[0]);
    multi = atoi(arg[1]);
    nthresh = atoi(arg[2]);
    if (pfraction < 0.0 || pfraction > 1.0) 
      error->all(FLERR,"Illegal pin command");
    if (multi != 0 && multi != 1) error->all(FLERR,"Illegal pin command");
    if (nthresh < 0) error->all(FLERR,"Illegal pin command");
    pin_create();
  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsStrainPin::grow_app()
{
  spin = iarray[0];
  strain = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsStrainPin::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
     // Allow for 'pinned' spin values of 'Q'
    if (spin[i] < 1 || spin[i] > nspins+1) flag = 1;
    if (strain[i] < 0.0) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppPottsStrainPin::site_propensity(int i)
{
  // no events for a pinned site
  if (spin[i] > nspins) return 0.0;

  // events = spin flips to neighboring site different than self
  // disallow flip to pinned site
  // disallow wild flips = flips to value different than all neighs

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == spin[i] || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  // for each flip:
  // compute energy difference between initial and final state
  // include strain scaling on effective temperature
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double efinal,scale;
  double prob = 0.0;

  for (m = 0; m < nevent; m++) {
    spin[i] = unique[m];
    efinal = site_energy(i);
    if (efinal <= einitial) {
      scale = 1.0 + strain[i];
      prob += 1.0/scale;
    } else if (temperature > 0.0) {
      scale = 1.0 + strain[i];
      prob += exp((einitial-efinal)*(t_inverse/scale));
    }
  }

  spin[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPottsStrainPin::site_event(int i, RandomPark *random)
{
  int j,m,value;

  // pick one event from total propensity by accumulating its probability
  // disallow flip to pinned site
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
    // Prevent 'flip' to 'pinned' neighbor or neighbor with same spin
    if (value == oldstate || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;

    spin[i] = value;
    efinal = site_energy(i);
    if (efinal <= einitial) {
      scale = 1.0 + strain[i];
      prob += 1.0/scale;
    } else if (temperature > 0.0) {
      scale = 1.0 + strain[i];
      prob += exp((einitial-efinal)*(t_inverse/scale));
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

void *AppPottsStrainPin::extract_app(char *name)
{
  if (strcmp(name,"nspins") == 0) return (void *) &nspins;
  if (strcmp(name,"strain") == 0) return (void *) strain;
  return NULL;
}

/* ----------------------------------------------------------------------
   change some sites to pinned sites
   user params = pfraction, multi, nthresh
------------------------------------------------------------------------- */

void AppPottsStrainPin::pin_create()
{
  int i,j,m,nattempt,nme,npin,ndiff;
  int flag,flags[2],flagall[2];
  tagint iglobal;

  int ndesired = static_cast<int> (pfraction*nglobal);

  RandomPark *random = new RandomPark(ranmaster->uniform());

  // single site inclusions
  // only put local sites into hash
  // nthresh = 0 for insertion anywhere
  // nthresh > 0 for insertion at grain boundaries

  if (!multi) {
    std::map<tagint,int> hash;
    for (i = 0; i < nlocal; i++)
      hash.insert(std::pair<tagint,int> (id[i],i));
    std::map<tagint,int>::iterator loc;

    npin = 0;
    while (npin < ndesired) {
      nattempt = ndesired - npin;
      for (i = 0; i < nattempt; i++) {
	iglobal = random->tagrandom(nglobal);
	loc = hash.find(iglobal);
	if (loc != hash.end()) {
	  if (nthresh == 0) spin[loc->second] = nspins+1;
	  else {
	    m = loc->second;
	    ndiff = 0;
	    for (j = 0; j < numneigh[m]; j++)
	      if (spin[m] != spin[neighbor[m][j]]) ndiff++;
	    if (ndiff >= nthresh) spin[m] = nspins+1;
	  }
	}
      }
      
      nme = 0;
      for (i = 0; i < nlocal; i++)
	if (spin[i] > nspins) nme++;
      MPI_Allreduce(&nme,&npin,1,MPI_INT,MPI_SUM,world);
    }

  // multi site inclusions
  // put local and ghost sites into hash
  // nthresh = 0 for insertion anywhere
  // nthresh > 0 for insertion at grain boundaries

  } else if (multi) {
    std::map<tagint,int> hash;
    for (i = 0; i < nlocal+nghost; i++)
      hash.insert(std::pair<tagint,int> (id[i],i));
    std::map<tagint,int>::iterator loc;

    tagint *list = new tagint[maxneigh+1];

    npin = 0;
    while (npin < ndesired) {
      iglobal = random->tagrandom(nglobal);
      loc = hash.find(iglobal);
      if (loc != hash.end() && loc->second < nlocal) {
	flag = 1;
	i = loc->second;
	if (spin[i] > nspins) flag = 0;
	for (j = 0; j < numneigh[i]; j++)
	  if (spin[neighbor[i][j]] > nspins) flag = 0;
	if (nthresh) {
	  ndiff = 0;
	  for (j = 0; j < numneigh[i]; j++)
	    if (spin[i] != spin[neighbor[i][j]]) ndiff++;
	  if (ndiff < nthresh) flag = 0;
	}

	if (flag) {
	  flags[0] = me+1;
	  flags[1] = numneigh[i] + 1;
	  spin[i] = nspins+1;
	  for (j = 0; j < numneigh[i]; j++) {
	    spin[neighbor[i][j]] = nspins+1;
	    list[j] = id[neighbor[i][j]];
	  }
	  list[j++] = id[i];
	} else flags[0] = flags[1] = 0;
      } else flags[0] = flags[1] = 0;
      
      MPI_Allreduce(&flags,&flagall,2,MPI_INT,MPI_SUM,world);

      if (flagall[0]) {
	MPI_Bcast(list,flagall[1],MPI_INT,flagall[0]-1,world);
	for (i = 0; i < flagall[1]; i++) {
	  loc = hash.find(list[i]);
	  if (loc != hash.end()) spin[loc->second] = nspins+1;
	}

	nme = 0;
	for (i = 0; i < nlocal; i++)
	  if (spin[i] > nspins) nme++;
	MPI_Allreduce(&nme,&npin,1,MPI_INT,MPI_SUM,world);
      }
    }

    delete [] list;
  }

  delete random;
}

/* ----------------------------------------------------------------------
   push new site onto stack and assign new id
 ------------------------------------------------------------------------- */

void AppPottsStrainPin::push_new_site(int i, int* cluster_ids, int id,
                                      std::stack<int>* cluststack)
{
  int isite = spin[i];

  if (isite != nspins+1) {
    cluststack->push(i);
    cluster_ids[i] = id;
  }
}
