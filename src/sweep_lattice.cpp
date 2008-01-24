/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sweep_lattice.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS;

#define DELTA 100

/* ---------------------------------------------------------------------- */

SweepLattice::SweepLattice(SPK *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  
  // parse optional args

  Lmask = false;
  mask = NULL;
  Lstrict = false;
  Lpicklocal = false;
  Lkmc = false;
  ranlat = NULL;
  delt = 1.0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"strict") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lstrict = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lstrict = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mask") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lmask = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lmask = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"picklocal") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lpicklocal = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lpicklocal = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kmc") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Lkmc = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Lkmc = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      delt = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal sweep_style command");
  }
  // initializations

  int nsectormax = 8;
  for (int i = 0; i < nsectormax; i++) {
    sector[i].solve = NULL;
    sector[i].propensity = NULL;
    sector[i].i2site = NULL;
    sector[i].site2i = NULL;
    sector[i].sites = NULL;
    sector[i].border = NULL;
  }

  // communicator needed between sweep sectors

  comm = new CommLattice(spk);
}

/* ---------------------------------------------------------------------- */

SweepLattice::~SweepLattice()
{
  delete random;
  delete comm;

  memory->sfree(mask);
  memory->sfree(ranlat);

  int nsectormax = 8;
  for (int isector = 0; isector < nsectormax; isector++) {
    delete sector[isector].solve;
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].i2site);
    memory->sfree(sector[isector].site2i);
    memory->sfree(sector[isector].sites);
    memory->sfree(sector[isector].border);
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice::init()
{
  int i,j,k,m;

  applattice = (AppLattice *) app;

  int sitecustom = applattice->sitecustom;
  lattice = applattice->lattice;

  dimension = applattice->dimension;
  nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  int *id = applattice->id;
  double **xyz = applattice->xyz;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  temperature = applattice->temperature;
  if (temperature != 0.0) t_inverse = 1.0/temperature;

  // app-specific settings

  masklimit = applattice->masklimit;
  delghost = applattice->delghost;
  dellocal = applattice->dellocal;

  delcol = delghost+1;
  if (Lstrict) ncolor = delcol*delcol;
  else ncolor = 1;

  // setup sectors
  // loop over all ownded sites to assign each to a sector (4/8 in 2d/3d)
  // create site2i array for each sector

  if (dimension == 2) nsector = 4;
  else if (dimension == 3) nsector = 8;

  for (int isector = 0; isector < nsector; isector++) {
    sector[isector].nlocal = 0;
    sector[isector].nmax = 0;
    delete sector[isector].solve;
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].i2site);
    memory->sfree(sector[isector].site2i);
    memory->sfree(sector[isector].sites);
    memory->sfree(sector[isector].border);
    sector[isector].solve = NULL;
    sector[isector].propensity = NULL;
    sector[isector].i2site = NULL;
    sector[isector].site2i = NULL;
    sector[isector].sites = NULL;
    sector[isector].border = NULL;
  }

  int iwhich,jwhich,kwhich,isector;
  double xmid = 0.5 * (applattice->subxlo + applattice->subxhi);
  double ymid = 0.5 * (applattice->subylo + applattice->subyhi);
  double zmid = 0.5 * (applattice->subzlo + applattice->subzhi);

  for (i = 0; i < nlocal; i++) {
    if (xyz[i][0] < xmid) iwhich = 0;
    else iwhich = 1;
    if (xyz[i][1] < ymid) jwhich = 0;
    else jwhich = 1;
    if (xyz[i][2] < zmid) kwhich = 0;
    else kwhich = 1;

    if (dimension == 2) isector = 2*iwhich + jwhich;
    else isector = 4*iwhich + 2*jwhich + kwhich;

    if (sector[isector].nlocal == sector[isector].nmax) {
      sector[isector].nmax += DELTA;
      sector[isector].site2i = 
	(int *) memory->srealloc(sector[isector].site2i,
				 sector[isector].nmax*sizeof(int),
				 "applattice:site2i");
    }
    sector[isector].site2i[sector[isector].nlocal++] = i;
  }

  // init communication for ghost sites

  comm->init(this,delghost,dellocal,NULL);

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    mask = (char *) memory->smalloc(nlocal*sizeof(char),"sweeplattice:mask");
    for (i = 0; i < nlocal; i++) mask[i] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    ranlat = (RandomPark *) memory->smalloc(nlocal*sizeof(RandomPark *),
					    "sweeplattice:ranlat");
    for (i = 0; i < nlocal; i++)
      ranlat[i].init(seed+id[i]);
  }

  // setup KMC solver and propensity arrays for each sector
  // i2site = mapping to sector site for i owned by sector, else -1
  // nborder = # of sites in sector with a neighbor who is outside sector
  // border = list of these sites, stored as lattice indices
  // propensity init requires ghost cell info for entire sub-domain

  if (Lkmc) {
    comm->all();

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nlocal;
      sector[isector].propensity = 
	(double *) memory->smalloc(nsites*sizeof(double),"sweep:propensity");

      sector[isector].i2site = 
	(int *) memory->smalloc(ntotal*sizeof(int),"sweep:i2site");
      for (i = 0; i < ntotal; i++) sector[isector].i2site[i] = -1;
      for (m = 0; m < nsites; m++) 
	sector[isector].i2site[sector[isector].site2i[m]] = m;

      sector[isector].nborder = 
	find_border_sites(nsites,sector[isector].site2i,
			  numneigh,neighbor,&sector[isector].border);
      sector[isector].sites =
	(int *) memory->smalloc(sector[isector].nborder*sizeof(int),
			       "sweep:sites");

      for (m = 0; m < nsites; m++)
	sector[isector].propensity[m] =
	  applattice->site_propensity(sector[isector].site2i[m],0);

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }

  // set sweep method function ptr

  if (Lkmc) {
    if (Lmask || Lpicklocal)
      error->all("Combination of sweep flags is unsupported");
    sweeper = &SweepLattice::sweep_sector_kmc;
  } else if (Lstrict) {
    if (Lmask) {
      if (Lpicklocal || sitecustom)
	error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice::sweep_sector_mask_strict;
    } else {
      if (Lpicklocal || sitecustom)
	error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice::sweep_sector_strict;
    }
  } else {
    if (Lmask) {
      if (sitecustom)
	error->all("Combination of sweep flags is unsupported");
      else if (Lpicklocal)
	sweeper = &SweepLattice::sweep_sector_mask_picklocal;
      else sweeper = &SweepLattice::sweep_sector_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else if (sitecustom == 0)
	sweeper = &SweepLattice::sweep_sector_lattice;
      else
	sweeper = &SweepLattice::sweep_sector_data;
    }
  }
}

/* ----------------------------------------------------------------------
   create list of border sites for a sector
   border site = owned site in sector which has a neighbor outside the sector
 ------------------------------------------------------------------------- */

int SweepLattice::find_border_sites(int nsites, int *site2i,
				    int *numneigh, int **neighbor,
				    int **pborder)
{
  int i,j,m;

  int *border = *pborder;
  int nborder = 0;
  int nmax = 0;

  // put lattice index of each site into hash map
  // use lattice index since will be looking up neighbor[][] which stores
  //   a lattice index

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;
  for (int m = 0; m < nsites; m++)
    hash.insert(std::pair<int,int> (site2i[m],0));
  
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    for (j = 0; j < numneigh[i]; j++) {
      if (hash.find(neighbor[i][j]) == hash.end()) break;
    }
    if (j == numneigh[i]) continue;
    
    if (nborder == nmax) {
      nmax += DELTA;
      border = (int *) memory->srealloc(border,nmax*sizeof(int),
					"applattice:border");
    }
    border[nborder++] = i;
  }

  *pborder = border;
  return nborder;
}
  
/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice::do_sweep(double &dt)
{
  for (int icolor = 0; icolor < ncolor; icolor++)
    for (int isector = 0; isector < nsector; isector++) {
      timer->stamp();
      comm->sector(isector);
      timer->stamp(TIME_COMM);

      (this->*sweeper)(icolor,isector);

      timer->stamp(TIME_SOLVE);

      comm->reverse_sector(isector);
      timer->stamp(TIME_COMM);
    }

  dt = delt;
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites using lattice array as state
   application picks a new state for the site
   compute energy change due to state change
   no energy change = accept
   downhill energy change = accept
   uphill energy change = accept/reject via Boltzmann factor
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_lattice(int icolor, int isector)
{
  int i,j,m,oldstate;
  double einitial,efinal;

  int *site2i = sector[isector].site2i;
  int nlocal = sector[isector].nlocal;

  for (m = 0; m < nlocal; m++) {
    i = site2i[m];
    oldstate = lattice[i];
    einitial = applattice->site_energy(i);
  
    applattice->site_pick_random(i,random->uniform());
    efinal = applattice->site_energy(i);
    
    if (efinal <= einitial) continue;
    else if (temperature == 0.0) lattice[i] = oldstate;
    else if (random->uniform() > exp((einitial-efinal)*t_inverse))
      lattice[i] = oldstate;
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites using general data for state
   application picks a new state for the site
   compute energy change due to state change
   no energy change = accept
   downhill energy change = accept
   uphill energy change = accept/reject via Boltzmann factor
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_data(int icolor, int isector)
{
  int i,j,m;
  double einitial,efinal;

  int *site2i = sector[isector].site2i;
  int nlocal = sector[isector].nlocal;

  for (m = 0; m < nlocal; m++) {
    i = site2i[m];
    einitial = applattice->site_energy(i);
    applattice->site_save(i);
  
    applattice->site_pick_random(i,random->uniform());
    efinal = applattice->site_energy(i);
    
    if (efinal <= einitial) continue;
    else if (temperature == 0.0) applattice->site_restore(i);
    else if (random->uniform() > exp((einitial-efinal)*t_inverse))
      applattice->site_restore(i);
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   skip sites that can't change via mask
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_mask(int icolor, int isector)
{
}

/* ---------------------------------------------------------------------- */

void SweepLattice::sweep_sector_mask_picklocal(int icolor, int isector)
{
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_strict(int icolor, int isector)
{
}

/* ---------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_mask_strict(int icolor, int isector)
{
}

/* ----------------------------------------------------------------------
   generate events in each sector using KMC solver
 ------------------------------------------------------------------------- */

void SweepLattice::sweep_sector_kmc(int icolor, int isector)
{
  double dt,time;
  int done,isite,i;

  timer->stamp();

  // extract sector specific info from struct

  int *border = sector[isector].border;
  int nborder = sector[isector].nborder;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int *i2site = sector[isector].i2site;
  int *site2i = sector[isector].site2i;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // i2site, solver, propensity array

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  double *hold_propensity = applattice->propensity;
  applattice->propensity = propensity;
  int *hold_i2site = applattice->i2site;
  applattice->i2site = i2site;

  // update owned propensities for sites which neighbor a sector ghost
  // necessary since sector ghosts may have changed

  int nsites = 0;

  for (int m = 0; m < nborder; m++) {
    i = border[m];
    isite = i2site[i];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i,0);
  }

  solve->update(nsites,sites,propensity);

  // attribute this chunk of time to Communication,
  // because it is due to the spatial decomposition
  timer->stamp(TIME_COMM);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    time += dt;	
    if (isite < 0 || time >= delt) done = 1;
    else {
      i = site2i[isite];
      applattice->site_event(i,0);
      timer->stamp(TIME_APP);
    }
  }
 
  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = hold_propensity;
  applattice->i2site = hold_i2site;
}

/* ----------------------------------------------------------------------
   unset all masks on boundary
   may be out of date, due to state change on neighboring processor
   could reverse comm mask values, but that might be slower
 ------------------------------------------------------------------------- */

void SweepLattice::boundary_clear_mask()
{
}
