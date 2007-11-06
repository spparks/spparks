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

using namespace SPPARKS;

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

  // set sweep method function ptr

  if (Lkmc) {
    if (Lmask || Lpicklocal)
      error->all("Combination of sweep flags is unsupported");
    sweeper = &SweepLattice::sweep_sector_kmc;
  } else if (Lstrict) {
    if (Lmask) {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice::sweep_sector_mask_strict;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice::sweep_sector_strict;
    }
  } else {
    if (Lmask) {
      if (Lpicklocal) sweeper = &SweepLattice::sweep_sector_mask_picklocal;
      else sweeper = &SweepLattice::sweep_sector_mask;
    } else {
      if (Lpicklocal) error->all("Combination of sweep flags is unsupported");
      else sweeper = &SweepLattice::sweep_sector;
    }
  }

  // initializations

  int nsectormax = 8;
  for (int i = 0; i < nsectormax; i++) {
    sector[i].solve = NULL;
    sector[i].propensity = NULL;
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

  if (Lkmc) {
    int nsectormax = 8;
    for (int i = 0; i < nsectormax; i++) {
      delete sector[i].solve;
      memory->sfree(sector[i].propensity);
    }
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice::init()
{
  applattice = (AppLattice *) app;

  lattice = applattice->lattice;
  int *id = applattice->id;

  nlocal = applattice->nlocal;
  nghost = applattice->nghost;

  dimension = applattice->dimension;
  int procwest = applattice->procwest;
  int proceast = applattice->proceast;
  int procsouth = applattice->procsouth;
  int procnorth = applattice->procnorth;
  int procdown,procup;
  if (dimension == 3) {
    procdown = applattice->procdown;
    procup = applattice->procup;
  }

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

  if (dimension == 2) nsector = 4;
  else if (dimension == 3) nsector = 8;

  for (int isector = 0; isector < nsector; isector++) {
    sector[isector].n = 1;
    delete sector[isector].solve;
    memory->sfree(sector[isector].propensity);
  }

  // init communication for ghost sites

  if (dimension == 2)
    comm->init(nlocal,procwest,proceast,procsouth,procnorth,
	       delghost,dellocal);
  else if (dimension == 3)
    comm->init(nlocal,procwest,proceast,procsouth,procnorth,procdown,procup,
	       delghost,dellocal);

  // setup mask array
  // owned and ghost values referenced in app::site_clear_mask()

  if (Lmask && mask == NULL) {
    mask = (char *) memory->smalloc(nlocal*sizeof(char),"sweeplattice:mask");
    for (int i = 0; i < nlocal; i++) mask[i] = 0;
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    ranlat = (RandomPark *) memory->smalloc(nlocal*sizeof(RandomPark *),
					    "sweeplattice:ranlat");
    int isite;
    for (int i = 0; i < nlocal; i++)
      ranlat[i].init(seed+id[i]);
  }

  // setup KMC solver and propensity arrays for each sector
  // propensity init requires ghost cell info for entire sub-domain

  if (Lkmc) {
    int i,j,m;

    comm->all(lattice);

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = 1;
      sector[isector].propensity = 
	(double*) memory->smalloc(nsites*sizeof(double),"sweep:propensity");

      // init propensity values for this sector

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice::do_sweep(double &dt)
{
  dt = delt;
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   application picks a new state for the site
   compute energy change due to state change
   no energy change = accept
   downhill energy change = accept
   uphill energy change = accept/reject via Boltzmann factor
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector(int icolor, int isector)
{
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
}

/* ----------------------------------------------------------------------
   unset all masks on boundary
   may be out of date, due to state change on neighboring processor
   could reverse comm mask values, but that might be slower
 ------------------------------------------------------------------------- */

void SweepLattice::boundary_clear_mask()
{
}
