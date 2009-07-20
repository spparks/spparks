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
#include "app_diffusion2.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED,TOP};
enum{NO_ENERGY,LINEAR,LINEAR_TABLE,NONLINEAR};
enum{NO_GEOMETRY,SURFACE,VOID_FRACTION,PORE};
enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};  // from app_lattice.cpp

#define DELTAEVENT 100000

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppDiffusion2::AppDiffusion2(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;     // change to 3 for nonlinear, see below
  allow_rejection = 0;   // change to 1 for linear, see below
  allow_masking = 0;

  // parse arguments

  if (narg < 4) error->all("Illegal app_style command");

  if (strcmp(arg[1],"off") == 0) engflag = NO_ENERGY;
  else if (strcmp(arg[1],"linear") == 0) engflag = LINEAR;
  else if (strcmp(arg[1],"linear/table") == 0) engflag = LINEAR_TABLE;
  else if (strcmp(arg[1],"nonlinear") == 0) engflag = NONLINEAR;
  else error->all("Illegal app_style command");

  if (engflag == NONLINEAR) delpropensity = 3;
  if (engflag == LINEAR) allow_rejection = 1;

  int input_style,iarg;
  double level,fraction,xc,yc,zc,diameter,thickness;

  if (strcmp(arg[2],"none") == 0) {
    if (narg < 3) error->all("Illegal app_style command");
    input_style = NO_GEOMETRY;
    iarg = 3;
  } else if (strcmp(arg[2],"surface") == 0) {
    if (narg < 4) error->all("Illegal app_style command");
    input_style = SURFACE;
    level = atof(arg[3]);
    if (level < 0.0) error->all("Illegal app_style command");
    iarg = 4;
  } else if (strcmp(arg[2],"void") == 0) {
    if (narg < 4) error->all("Illegal app_style command");
    input_style = VOID_FRACTION;
    fraction = atof(arg[3]);
    if (fraction <= 0.0 || fraction >= 1.0)
      error->all("Illegal app_style command");
    iarg = 4;
  } else if (strcmp(arg[2],"pore") == 0) {
    if (narg < 8) error->all("Illegal app_style command");
    input_style = PORE;
    xc = atof(arg[3]);
    yc = atof(arg[4]);
    zc = atof(arg[5]);
    diameter = atof(arg[6]);
    thickness = atof(arg[7]);
    iarg = 8;
  } else error->all("Illegal app_style command");

  options(narg-iarg,&arg[iarg]);

  // define lattice and partition it across processors

  create_lattice();

  // for linear and linear/table:
  // make sites large enough for 2 sites and 1st/2nd nearest neighbors
  // for nonlinear:
  // make sites large enough for 2 sites and their 1,2,3 nearest neighbors

  if (engflag == LINEAR || engflag == LINEAR_TABLE) {
    esites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
    psites = NULL;
  } else if (engflag == NONLINEAR) {
    int nmax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
    esites = new int[2*nmax];
    psites = new int[2*nmax];
  } else esites = psites = NULL;

  echeck = pcheck = NULL;

  // event list

  events = NULL;
  firstevent = NULL;

  // default settings for app-specific commands

  depflag = 0;
  ndeposit = ndeposit_failed = 0;

  ecoord = new double[maxneigh+1];
  for (int i = 0; i <= maxneigh; i++) ecoord[i] = 0.0;

  barrierflag = 0;
  barrier = 
    memory->create_2d_double_array(maxneigh+1,maxneigh+1,"app:barrier");
  for (int i = 0; i <= maxneigh; i++)
    for (int j = 0; j <= maxneigh; j++)
      barrier[i][j] = 0.0;

  // sweeping timestep

  dt_sweep = 1.0/maxneigh;

  // initialize my portion of lattice as file or SURFACE or VOID_FRACTION
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  RandomPark *random = new RandomPark(ranmaster->uniform());

  if (infile) read_file();

  // SURFACE is sites < level = OCCUPIED, rest VACANT, except at TOP

  else if (input_style == SURFACE) {
    std::map<int,int> hash;
    for (int i = 0; i < nlocal; i++)
      hash.insert(std::pair<int,int> (id[i],i));
    std::map<int,int>::iterator loc;
    
    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      loc = hash.find(iglobal);
      if (loc != hash.end()) {
	if (dimension == 2) {
	  if (xyz[loc->second][1] <= level) {
	    lattice[loc->second] = OCCUPIED;
	  } else if (xyz[loc->second][1] >= yprd-2.0*latconst)
	    lattice[loc->second] = TOP;
	  else lattice[loc->second] = VACANT;
	} else {
	  if (xyz[loc->second][2] <= level)
	    lattice[loc->second] = OCCUPIED;
	  else if (xyz[loc->second][2] >= zprd-2.0*latconst)
	    lattice[loc->second] = TOP;
	  else lattice[loc->second] = VACANT;
	}
      }
    }

  // VOID_FRACTION is sites = VACANT/OCCUPIED with fraction OCCUPIED

  } else if (input_style == VOID_FRACTION) {
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

  // PORE is sites = VACANT/OCCUPIED
  // pore is in thin film at (xc,yc,zc) with diameter and thickness

  } else if (input_style == PORE) {
    double x,y,z;
    int isite;
    for (int i = 0; i < nlocal; i++) {
      x = xyz[i][0];
      y = xyz[i][1];
      z = xyz[i][2];
      if (dimension == 2) {
	if (y > yc + 0.5*thickness || y < yc - 0.5*thickness) isite = VACANT;
	else isite = OCCUPIED;
	if (isite == OCCUPIED) {
	  if ((x-xc)*(x-xc) < 0.25*diameter*diameter) isite = VACANT;
	}
      } else {
	if (z > zc + 0.5*thickness || z < zc - 0.5*thickness) isite = VACANT;
	else isite = OCCUPIED;
	if (isite == OCCUPIED) {
	  if ((x-xc)*(x-xc) + (y-yc)*(y-yc) < 0.25*diameter*diameter)
	    isite = VACANT;
	}
      }
      lattice[i] = isite;
    }
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

AppDiffusion2::~AppDiffusion2()
{
  delete [] esites;
  delete [] psites;
  delete [] echeck;
  delete [] pcheck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->destroy_2d_double_array(barrier);
  delete [] ecoord;
}

/* ---------------------------------------------------------------------- */

void AppDiffusion2::input_app(char *command, int narg, char **arg)
{
  double PI = 4.0*atan(1.0);

  if (strcmp(command,"ecoord") == 0) {
    if (engflag != NONLINEAR)
      error->all("Can only use ecoord command with "
		 "app_style diffusion nonlinear");
    if (narg != 2) error->all("Illegal ecoord command");

    int lo,hi;
    bounds(arg[0],maxneigh,lo,hi);
    double value = atof(arg[1]);

    for (int i = lo; i <= hi; i++) ecoord[i] = value;

  } else if (strcmp(command,"deposition") == 0) {
    if (narg < 1) error->all("Illegal deposition command");
    if (strcmp(arg[0],"off") == 0) {
      if (narg != 1 ) error->all("Illegal deposition command");
      depflag = 0;
      return;
    }

    if (narg != 7) error->all("Illegal deposition command");
    depflag = 1;
    deprate = atof(arg[0]);
    dir[0] = atof(arg[1]);
    dir[1] = atof(arg[2]);
    dir[2] = atof(arg[3]);
    d0 = atof(arg[4]);
    coordlo = atoi(arg[5]);
    coordhi = atoi(arg[6]);
    if (deprate < 0.0) error->all("Illegal deposition command");
    if (dimension == 2 && (dir[1] >= 0.0 || dir[2] != 0.0))
      error->all("Illegal deposition command");
    if (dimension == 3 && dir[2] >= 0.0)
      error->all("Illegal deposition command");
    if (d0 < 0.0) error->all("Illegal deposition command");
    if (coordlo < 0 || coordhi > maxneigh || coordlo > coordhi)
      error->all("Illegal deposition command");

    double len = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] /= len;
    dir[1] /= len;
    dir[2] /= len;

  } else if (strcmp(command,"barrier") == 0) {
    if (narg > 3) error->all("Illegal barrier command");
    barrierflag = 1;

    if (narg == 1) {
      double q = atof(arg[0]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
	for (j = 0; j <= maxneigh; j++)
	  barrier[i][j] = q;
      if (q == 0.0) barrierflag = 0;

    } else if (narg == 2) {
      int delta = atoi(arg[0]);
      double q = atof(arg[1]);
      int i,j;
      for (i = 0; i <= maxneigh; i++)
	for (j = 0; j <= maxneigh; j++)
	  if (j-i == delta) barrier[i][j] = q;

    } else {
      int ilo,ihi,jlo,jhi;
      bounds(arg[0],maxneigh,ilo,ihi);
      bounds(arg[1],maxneigh,jlo,jhi);
      double q = atof(arg[2]);

      for (int i = ilo; i <= ihi; i++)
	for (int j = jlo; j <= jhi; j++)
	  barrier[i][j] = q;
    }

  } else error->all("Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppDiffusion2::init_app()
{
  // rejection and deposition are incompatible
  // no parallel deposition for now

  if (depflag && sweepflag != NOSWEEP)
    error->all("Cannot run rejection KMC with deposition");

  if (depflag && nprocs > 1)
    error->all("Cannot run deposition in parallel for now");

  delete [] echeck;
  echeck = new int[nlocal+nghost];
  delete [] pcheck;
  pcheck = new int[nlocal+nghost];

  if (engflag == LINEAR_TABLE || engflag == NONLINEAR) {
    memory->sfree(events);
    memory->sfree(firstevent);

    events = NULL;
    maxevent = 0;
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
  }
}

/* ---------------------------------------------------------------------- */

void AppDiffusion2::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = pcheck[i] = 0;

  if (engflag == LINEAR_TABLE || engflag == NONLINEAR) {
    nevents = 0;
    for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppDiffusion2::site_energy(int i)
{
  // energy only non-zero for OCCUPIED sites when energy included in model

  if (lattice[i] != OCCUPIED || engflag == NO_ENERGY) return 0.0;

  // energy is a non-linear function of coordination number
  // calculate from user-specified tabulated values

  if (engflag == NONLINEAR) {
    int n = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) n++;
    return ecoord[n];
  }

  // energy is a linear function of coordination number, just count bonds

  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
   null bin extends to size maxneigh
------------------------------------------------------------------------- */

void AppDiffusion2::site_event_rejection(int i, RandomPark *random)
{
  double einitial,edelta;

  // OCCUPIED site exchanges with random neighbor if VACANT

  if (lattice[i] != OCCUPIED) return;
  int iran = (int) (maxneigh*random->uniform());
  if (iran > maxneigh) iran = maxneigh-1;
  int j = neighbor[i][iran];
  if (lattice[j] != VACANT) return;

  // accept or reject via energy and barrier model
  // factor of 2 in edelta accounts for energy change of neighbors of I,J

  int hop = 0;
  if (engflag != NO_ENERGY) einitial = site_energy(i);

  lattice[i] = VACANT;
  lattice[j] = OCCUPIED;

  if (engflag == NO_ENERGY) {
    if (!barrierflag) hop = 1;
    else if (temperature > 0.0) {
      if (random->uniform() < exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse))
	  hop = 1;
    }

  } else {
    edelta = site_energy(j) - einitial;

    if (!barrierflag) {
      if (edelta <= 0.0) hop = 1;
      else if (temperature > 0.0) {
	if (random->uniform() < exp(-2.0*edelta*t_inverse)) hop = 1;
      }
    } else if (temperature > 0.0) {
      if (edelta <= 0.0) {
	if (random->uniform() < 
	    exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse)) hop = 1;
      } else {
	if (random->uniform() < 
	    exp((-2.0*edelta-barrier[ncoord(i)-1][ncoord(j)]) * t_inverse))
	  hop = 1;
      }
    }
  }

  if (hop) naccept++;
  else {
    lattice[i] = OCCUPIED;
    lattice[j] = VACANT;
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppDiffusion2::site_propensity(int i)
{
  if (engflag == LINEAR) return site_propensity_linear(i);
  else if (engflag == LINEAR_TABLE) return site_propensity_table(i);
  else return site_propensity_nonlinear(i);
}

/* ---------------------------------------------------------------------- */

double AppDiffusion2::site_propensity_linear(int i)
{
  int j;
  double einitial,edelta,prob;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engflag, compute edelta between initial and final state
  // factor of 2 in edelta accounts for energy change of neighbors of I,J
  // if barrierflag, compute coordination of sites I,J, with hop atom removed
  // propensity is function of edelta, barrier, up/down hill, temperature

  if (lattice[i] != OCCUPIED) return 0.0;

  prob = 0.0;
  if (engflag != NO_ENERGY) einitial = site_energy(i);

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANT) {
      lattice[i] = VACANT;
      lattice[j] = OCCUPIED;

      if (engflag == NO_ENERGY) {
	if (!barrierflag) prob += 1.0;
	else if (temperature > 0.0)
	  prob += exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);

      } else {
	edelta = site_energy(j) - einitial;

	if (!barrierflag) {
	  if (edelta <= 0.0) prob += 1.0;
	  else if (temperature > 0.0) 
	    prob += exp(-2.0*edelta*t_inverse);
	} else if (temperature > 0.0) {
	  if (edelta <= 0.0)
	    prob += exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);
	  else
	    prob += exp((-2.0*edelta-barrier[ncoord(i)-1][ncoord(j)]) * 
			t_inverse);
	}
      }

      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) prob += deprate;

  return prob;
}

/* ---------------------------------------------------------------------- */

double AppDiffusion2::site_propensity_table(int i)
{
  int j;
  double einitial,edelta,probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engflag, compute edelta between initial and final state
  // factor of 2 in edelta accounts for energy change of neighbors of I,J
  // if barrierflag, compute coordination of sites I,J, with hop atom removed
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) return 0.0;

  proball = 0.0;
  if (engflag != NO_ENERGY) einitial = site_energy(i);

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANT) {
      lattice[i] = VACANT;
      lattice[j] = OCCUPIED;
      probone = 0.0;

      if (engflag == NO_ENERGY) {
	if (!barrierflag) probone = 1.0;
	else if (temperature > 0.0)
	  probone = exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);

      } else {
	edelta = site_energy(j) - einitial;

	if (!barrierflag) {
	  if (edelta <= 0.0) probone = 1.0;
	  else if (temperature > 0.0) 
	    probone = exp(-2.0*edelta*t_inverse);
	} else if (temperature > 0.0) {
	  if (edelta <= 0.0)
	    probone = exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);
	  else
	    probone = exp((-2.0*edelta-barrier[ncoord(i)-1][ncoord(j)]) * 
			  t_inverse);
	}
      }

      if (probone > 0.0) {
	add_event(i,j,probone);
	proball += probone;
      }

      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate);
    proball += deprate;
  }

  return proball;
}

/* ---------------------------------------------------------------------- */

double AppDiffusion2::site_propensity_nonlinear(int i)
{
  int j,k,m,nsites;
  double einitial,efinal,edelta,probone,proball;

  // events = OCCUPIED site exchanges with adjacent VACANT site
  // if engflag, compute edelta between initial and final state
  // since eng is nonlinear, this must include eng of neighbor sites of I,J
  // if barrierflag, compute coordination of sites I,J, with hop atom removed
  // propensity is function of edelta, barrier, up/down hill, temperature

  clear_events(i);

  if (lattice[i] != OCCUPIED) return 0.0;

  proball = 0.0;

  for (int ineigh = 0; ineigh < numneigh[i]; ineigh++) {
    j = neighbor[i][ineigh];
    if (lattice[j] == VACANT) {
      probone = 0.0;

      if (engflag == NO_ENERGY) {
	lattice[i] = VACANT;
	lattice[j] = OCCUPIED;
	if (!barrierflag) probone = 1.0;
	else if (temperature > 0.0)
	  probone = exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);

      } else {

	// einitial = i,j and their neighbors
	// use pcheck[] to avoid recomputing energy of same site

	einitial = site_energy(i) + site_energy(j);
	nsites = 0;
	psites[nsites++] = i;
	psites[nsites++] = j;
	pcheck[i] = pcheck[j] = 1;
	
	for (k = 0; k < numneigh[i]; k++) {
	  m = neighbor[i][k];
	  if (pcheck[m]) continue;
	  einitial += site_energy(m);
	  psites[nsites++] = m;
	  pcheck[m] = 1;
	}
	for (k = 0; k < numneigh[j]; k++) {
	  m = neighbor[j][k];
	  if (pcheck[m]) continue;
	  einitial += site_energy(m);
	  psites[nsites++] = m;
	  pcheck[m] = 1;
	}

	for (m = 0; m < nsites; m++) pcheck[psites[m]] = 0;

	lattice[i] = VACANT;
	lattice[j] = OCCUPIED;

	// efinal = i,j and their neighbors
	// use pcheck[] to avoid recomputing energy of same site

	efinal = site_energy(i) + site_energy(j);
	nsites = 0;
	psites[nsites++] = i;
	psites[nsites++] = j;
	pcheck[i] = pcheck[j] = 1;
	
	for (k = 0; k < numneigh[i]; k++) {
	  m = neighbor[i][k];
	  if (pcheck[m]) continue;
	  efinal += site_energy(m);
	  psites[nsites++] = m;
	  pcheck[m] = 1;
	}
	for (k = 0; k < numneigh[j]; k++) {
	  m = neighbor[j][k];
	  if (pcheck[m]) continue;
	  efinal += site_energy(m);
	  psites[nsites++] = m;
	  pcheck[m] = 1;
	}
	
	for (m = 0; m < nsites; m++) pcheck[psites[m]] = 0;

	edelta = efinal - einitial;

	if (!barrierflag) {
	  if (edelta <= 0.0) probone = 1.0;
	  else if (temperature > 0.0) 
	    probone = exp(-edelta*t_inverse);
	} else if (temperature > 0.0) {
	  if (edelta <= 0.0)
	    probone = exp(-barrier[ncoord(i)-1][ncoord(j)]*t_inverse);
	  else
	    probone = exp((-edelta-barrier[ncoord(i)-1][ncoord(j)]) * 
			  t_inverse);
	}
      }

      if (probone > 0.0) {
	add_event(i,j,probone);
	proball += probone;
      }
      
      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate);
    proball += deprate;
  }

  return proball;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppDiffusion2::site_event(int i, class RandomPark *random)
{
  if (engflag == LINEAR) return site_event_linear(i,random);
  else if (engflag == LINEAR_TABLE) return site_event_table(i,random);
  else return site_event_nonlinear(i,random);
}

/* ---------------------------------------------------------------------- */

void AppDiffusion2::site_event_linear(int i, class RandomPark *random)
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
      einitial = site_energy(i);
      lattice[i] = VACANT;
      lattice[j] = OCCUPIED;
      efinal = site_energy(j);
      if (efinal <= einitial) prob += 1.0;
      else if (temperature > 0.0) prob += exp(2.0*(einitial-efinal)*t_inverse);
      if (prob >= threshhold) break;
      lattice[i] = OCCUPIED;
      lattice[j] = VACANT;
    }
  }

  // deposition event
  // find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (depflag && i == 0 && prob < threshhold) {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = TOP;
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

/* ---------------------------------------------------------------------- */

void AppDiffusion2::site_event_table(int i, class RandomPark *random)
{
  int j,jj,k,kk,m,mm;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // if partner >= 0, hop event, else deposition event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  j = events[ievent].partner;
  if (j >= 0) {
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  } else {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = TOP;
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

/* ---------------------------------------------------------------------- */

void AppDiffusion2::site_event_nonlinear(int i, class RandomPark *random)
{
  int j,jj,k,kk,kkk,m,mm,mmm,isite;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // if partner >= 0, hop event, else deposition event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  j = events[ievent].partner;
  if (j >= 0) {
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;
  } else {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = TOP;
    lattice[m] = OCCUPIED;
    i = j = m;
  }

  if (site_energy(j) > 1500.0) exit(1);

  // compute propensity changes for self and swap site and their 1,2,3 neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
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
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
	  esites[nsites++] = isite;
	  echeck[isite] = 1;
	}
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
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
	  esites[nsites++] = isite;
	  echeck[isite] = 1;
	}
      }
    }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ---------------------------------------------------------------------- */

int AppDiffusion2::ncoord(int i)
{
  int count = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  return count;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppDiffusion2::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppDiffusion2::add_event(int i, int partner, double propensity)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].partner = partner;
  events[freeevent].next = firstevent[i];
  events[freeevent].propensity = propensity;
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom
------------------------------------------------------------------------- */

int AppDiffusion2::find_deposition_site(RandomPark *random)
{
  // pick a random position at top of box

  double start[3];
  start[0] = boxxlo + (boxxhi-boxxlo)*random->uniform();
  start[1] = boxylo + (boxyhi-boxylo)*random->uniform();
  if (dimension == 3) start[2] = boxzhi;
  else start[2] = 0.0;

  // for each vacant site:
  // discard site if neighbor count not between coordlo and coordhi
  // find site whose projected distance is closest to start point

  int i,j,ncount;
  double dist_projected;

  int closesite = -1;
  double closedist = 1.0e20;

  for (i = 0; i < nlocal; i++) {
    if (lattice[i] != VACANT) continue;
    ncount = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) ncount++;
    if (ncount < coordlo || ncount > coordhi) continue;
    if (exceed_limit(i,start,dist_projected)) continue;
    if (dist_projected < closedist) {
      closedist = dist_projected;
      closesite = i;
    }
  }

  if (closesite < 0) ndeposit_failed++;
  else ndeposit++;

  return closesite;
}

/* ----------------------------------------------------------------------
   compute projected_dist from site M to incident line of deposition
   return 1 if normal_dist exceeds d0 limit, else return 0
   must account for periodic images in XY of incident line and start point
------------------------------------------------------------------------- */

int AppDiffusion2::exceed_limit(int m, double *start, double &pdist)
{
  int increment,iprd,jprd;

  iprd = jprd = 0;
  double d0sq = d0*d0;

  double distsq = distsq_to_dir(m,start,iprd,jprd,pdist);
  double newdistsq = distsq_to_dir(m,start,iprd-1,jprd,pdist);
  if (newdistsq < distsq) increment = -1;
  else increment = 1;
  while ((newdistsq = distsq_to_dir(m,start,iprd,jprd,pdist)) < distsq) {
    distsq = newdistsq;
    iprd += increment;
  }
  iprd -= increment;

  if (dimension == 3) {
    newdistsq = distsq_to_dir(m,start,iprd,jprd-1,pdist);
    if (newdistsq < distsq) increment = -1;
    else increment = 1;
    while ((newdistsq = distsq_to_dir(m,start,iprd,jprd,pdist)) < distsq) {
      distsq = newdistsq;
      jprd += increment;
    }
  }
  jprd -= increment;

  if (distsq > d0sq) return 1;
  distsq = distsq_to_dir(m,start,iprd,jprd,pdist);
  return 0;
}

/* ----------------------------------------------------------------------
   compute projected_dist from site M to incident line of deposition
   return 1 if normal_dist exceeds d0 limit, else return 0
   must account for periodic images in XY of incident line and start point
------------------------------------------------------------------------- */

double AppDiffusion2::distsq_to_dir(int m, double *start,
				   int iprd, int jprd, double &pdist)
{
  double dot,distsq;
  double delta[3],projection[3],offset[3];

  delta[0] = xyz[m][0] + iprd*xprd - start[0];
  delta[1] = xyz[m][1] + jprd*yprd - start[1];
  delta[2] = xyz[m][2] - start[2];
    
  pdist = dir[0]*delta[0] + dir[1]*delta[1] + dir[2]*delta[2];
  projection[0] = pdist*dir[0];
  projection[1] = pdist*dir[1];
  projection[2] = pdist*dir[2];
  
  offset[0] = delta[0] - projection[0];
  offset[1] = delta[1] - projection[1];
  offset[2] = delta[2] - projection[2];
  return offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 0 to nmax,
     (3) i* = 0 to nmax, (4) *j = 0 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void AppDiffusion2::bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),0);
    nhi = MIN(atoi(str),nmax);
  } else if (strlen(str) == 1) {
    nlo = 0;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 0;
    nhi = MIN(atoi(ptr+1),nmax);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),0);
    nhi = nmax;
  } else {
    nlo = MAX(atoi(str),0);
    nhi = MIN(atoi(ptr+1),nmax);
  }
}
