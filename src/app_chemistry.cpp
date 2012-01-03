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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_chemistry.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "output.h"

using namespace SPPARKS_NS;

#define MAX_PRODUCT 5
#define AVOGADRO 6.023e23

/* ---------------------------------------------------------------------- */

AppChemistry::AppChemistry(SPPARKS *spk, int narg, char **arg) : 
  App(spk,narg,arg)
{
  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  // default settings

  volume = 0.0;

  nspecies = 0;
  sname = NULL;

  nreactions = 0;
  rname = NULL;
  nreactant = NULL;
  reactants = NULL;
  nproduct = NULL;
  products = NULL;
  rate = NULL;

  pcount = NULL;
  ndepends = NULL;
  depends = NULL;
  propensity = NULL;
  rcount = NULL;

  nevents = 0;
}

/* ---------------------------------------------------------------------- */

AppChemistry::~AppChemistry()
{
  for (int i = 0; i < nspecies; i++) delete [] sname[i];
  memory->sfree(sname);

  for (int i = 0; i < nreactions; i++) delete [] rname[i];
  memory->sfree(rname);

  memory->destroy(nreactant);
  memory->destroy(reactants);
  memory->destroy(nproduct);
  memory->destroy(products);
  memory->destroy(rate);

  memory->destroy(pcount);
  memory->destroy(ndepends);
  memory->destroy(depends);
  memory->destroy(propensity);
  memory->destroy(rcount);
}

/* ---------------------------------------------------------------------- */

void AppChemistry::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"add_reaction") == 0) add_reaction(narg,arg);
  else if (strcmp(command,"add_species") == 0) add_species(narg,arg);
  else if (strcmp(command,"count") == 0) set_count(narg,arg);
  else if (strcmp(command,"volume") == 0) set_volume(narg,arg);
  else error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppChemistry::init()
{
  // error check

  if (solve == NULL) error->all(FLERR,"No solver class defined");
  if (volume <= 0.0) error->all(FLERR,"Invalid volume setting");
  if (nreactions == 0)
    error->all(FLERR,"No reactions defined for chemistry app");

  factor_zero = AVOGADRO * volume;
  factor_dual = 1.0 / (AVOGADRO * volume);

  // determine reaction dependencies

  delete [] ndepends;
  memory->destroy(depends);
  memory->create(ndepends,nreactions,"chemistry:ndepends");
  build_dependency_graph();

  // zero reaction counts

  delete [] rcount;
  memory->create(rcount,nreactions,"chemistry:rcount");
  for (int m = 0; m < nreactions; m++) rcount[m] = 0;

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppChemistry::setup()
{
  // compute initial propensity for each reaction

  delete [] propensity;
  memory->create(propensity,nreactions,"chemistry:propensity");
  for (int m = 0; m < nreactions; m++) propensity[m] = compute_propensity(m);

  // initialize solver

  solve->init(nreactions,propensity);

  // setup of output

  nextoutput = output->setup(time);
}

/* ----------------------------------------------------------------------
   iterate on Gillespie solver
------------------------------------------------------------------------- */

void AppChemistry::iterate()
{
  int m,ireaction;
  double dt;

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    timer->stamp();
    ireaction = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    // check if solver failed to pick an event

    if (ireaction < 0) done = 1;
    else {

      // update particle counts due to reaction

      rcount[ireaction]++;
      for (m = 0; m < nreactant[ireaction]; m++)
	pcount[reactants[ireaction][m]]--;
      for (m = 0; m < nproduct[ireaction]; m++)
	pcount[products[ireaction][m]]++;

      // update propensities of dependent reactions
      // inform Gillespie solver of changes

      for (m = 0; m < ndepends[ireaction]; m++)
	propensity[depends[ireaction][m]] = 
	  compute_propensity(depends[ireaction][m]);
      solve->update(ndepends[ireaction],depends[ireaction],propensity);

      // update time by Gillepsie dt

      nevents++;
      time += dt;
      if (time >= stoptime) done = 1;
    }

    timer->stamp(TIME_APP);

    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppChemistry::stats(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10g %10d",time,nevents);
  strpnt += strlen(strpnt);

  for (int m = 0; m < nspecies; m++) {
    sprintf(strpnt," %d",pcount[m]);
    strpnt += strlen(strpnt);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppChemistry::stats_header(char *strtmp)
{
  char *strpnt = strtmp;
  sprintf(strpnt," %10s %10s","Time","Step");
  strpnt += strlen(strpnt);

  for (int m = 0; m < nspecies; m++) {
    sprintf(strpnt," %s",sname[m]);
    strpnt += strlen(strpnt);
  }
}

/* ---------------------------------------------------------------------- */

void AppChemistry::set_count(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal count command");
  
  int ispecies = find_species(arg[0]);
  if (ispecies < 0) {
    char *str = new char[128];
    sprintf(str,"Species ID %s does not exist",arg[0]);
    error->all(FLERR,str);
  }
  pcount[ispecies] = atoi(arg[1]);
}

/* ---------------------------------------------------------------------- */

void AppChemistry::add_reaction(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal reaction command");

  // store ID

  if (find_reaction(arg[0]) >= 0) {
    char *str = new char[128];
    sprintf(str,"Reaction ID %s already exists",arg[0]);
    error->all(FLERR,str);
  }

  int n = nreactions + 1;
  rname = (char **) memory->srealloc(rname,n*sizeof(char *),
					  "chemistry:rname");
  int nlen = strlen(arg[0]) + 1;
  rname[nreactions] = new char[nlen];
  strcpy(rname[nreactions],arg[0]);

  // grow reaction arrays

  memory->grow(nreactant,n,"chemistry:nreactnant");
  memory->grow(reactants,n,2,"chemistry:reactants");
  memory->grow(nproduct,n,"chemistry:nproduct");
  memory->grow(products,n,MAX_PRODUCT,"chemistry:products");
  memory->grow(rate,n,"chemistry:rate");

  // find which arg is numeric reaction rate

  char c;
  int iarg = 1;
  while (iarg < narg) {
    c = arg[iarg][0];
    if ((c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.') break;
    iarg++;
  }

  // error checks

  if (iarg == narg) error->all(FLERR,"Reaction has no numeric rate");
  if (iarg < 1 || iarg > 3) 
    error->all(FLERR,"Reaction must have 0,1,2 reactants");
  if (narg-1 - iarg > MAX_PRODUCT) 
    error->all(FLERR,"Reaction cannot have more than MAX_PRODUCT products");

  // extract reactant and product species names
  // if any species does not exist, create it

  nreactant[nreactions] = 0;
  for (int i = 1; i < iarg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    reactants[nreactions][i-1] = ispecies;
    nreactant[nreactions]++;
  }

  rate[nreactions] = atof(arg[iarg]);

  nproduct[nreactions] = 0;
  for (int i = iarg+1; i < narg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    products[nreactions][i - (iarg+1)] = ispecies;
    nproduct[nreactions]++;
  }
  
  nreactions++;
}

/* ---------------------------------------------------------------------- */

void AppChemistry::add_species(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal species command");

  // grow species arrays

  int n = nspecies + narg;
  sname = (char **) memory->srealloc(sname,n*sizeof(char *),
					  "chemistry:sname");
  pcount = (int *) memory->srealloc(pcount,n*sizeof(int),
					 "chemistry:pcount");

  for (int iarg = 0; iarg < narg; iarg++) {
    if (find_species(arg[iarg]) >= 0) {
      char *str = new char[128];
      sprintf(str,"Species ID %s already exists",arg[iarg]);
      error->all(FLERR,str);
    }
    int nlen = strlen(arg[iarg]) + 1;
    sname[nspecies+iarg] = new char[nlen];
    strcpy(sname[nspecies+iarg],arg[iarg]);
    pcount[nspecies+iarg] = 0;
  }
  nspecies += narg;
}

/* ---------------------------------------------------------------------- */

void AppChemistry::set_volume(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal volume command");
  volume = atof(arg[0]);
}

/* ----------------------------------------------------------------------
   return reaction index (0 to N-1) for a reaction ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppChemistry::find_reaction(char *str)
{
  for (int i = 0; i < nreactions; i++)
    if (strcmp(str,rname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return species index (0 to N-1) for a species ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppChemistry::find_species(char *str)
{
  for (int i = 0; i < nspecies; i++)
    if (strcmp(str,sname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   build dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void AppChemistry::build_dependency_graph()
{
  int i,j,k,m,n,mspecies,nspecies;

  // count the dependencies in flag array:
  // loop over reactants & products of each reaction
  // for each reaction m, mspecies = its reactants and products
  // for each species, loop over reactants of all reactions n
  // if a match, then set flag[n] since n is in dependency list of m

  int *flag = new int[nreactions];

  for (m = 0; m < nreactions; m++) {
    for (n = 0; n < nreactions; n++) flag[n] = 0;

    for (i = 0; i < nreactant[m]; i++) {
      mspecies = reactants[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) flag[n] = 1;
	}
      }
    }

    for (i = 0; i < nproduct[m]; i++) {
      mspecies = products[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) flag[n] = 1;
	}
      }
    }

    ndepends[m] = 0;
    for (n = 0; n < nreactions; n++) if (flag[n]) ndepends[m]++;
  }

  delete [] flag;

  // allocate depends array, 2nd dim is max of ndepends[]

  memory->destroy(depends);
  int nmax = 0;
  for (m = 0; m < nreactions; m++) nmax = MAX(nmax,ndepends[m]);
  memory->create(depends,nreactions,nmax,"chemistry:depends");

  // zero the dependencies

  for (m = 0; m < nreactions; m++) ndepends[m] = 0;

  // store the dependencies via same loops as before
  // k loop insures dependency was not already stored

  for (m = 0; m < nreactions; m++) ndepends[m] = 0;

  for (m = 0; m < nreactions; m++) {
    for (i = 0; i < nreactant[m]; i++) {
      mspecies = reactants[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) {
	    for (k = 0; k < ndepends[m]; k++)
	      if (n == depends[m][k]) break;
	    if (k == ndepends[m]) depends[m][ndepends[m]++] = n;
	  }
	}
      }
    }

    for (i = 0; i < nproduct[m]; i++) {
      mspecies = products[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) {
	    for (k = 0; k < ndepends[m]; k++)
	      if (n == depends[m][k]) break;
	    if (k == ndepends[m]) depends[m][ndepends[m]++] = n;
	  }
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute propensity of a single reaction
   for mono reaction: propensity = count * rate
   for dual reaction: propensity = count1 * count2 * rate / (Avogadro*Volume)
   for dual reaction: cut in half if reactants are same species
------------------------------------------------------------------------- */

double AppChemistry::compute_propensity(int m)
{
  double p;
  if (nreactant[m] == 0) p = factor_zero * rate[m];
  else if (nreactant[m] == 1) p = pcount[reactants[m][0]] * rate[m];
  else {
    if (reactants[m][0] == reactants[m][1]) 
      p = 0.5 * factor_dual * pcount[reactants[m][0]] * 
	(pcount[reactants[m][1]] - 1) * rate[m];
    else
      p = factor_dual * pcount[reactants[m][0]] * 
	pcount[reactants[m][1]] * rate[m];
  }
  return p;
}
