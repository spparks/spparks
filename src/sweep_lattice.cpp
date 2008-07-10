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

using namespace SPPARKS_NS;

#define DELTA 100

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SweepLattice::SweepLattice(SPPARKS *spk, int narg, char **arg) : 
  Sweep(spk,narg,arg)
{
  if (narg < 2) error->all("Illegal sweep_style command");

  seed = atoi(arg[1]);
  random = new RandomPark(seed);
  
  // parse optional args

  Lmask = false;
  mask = NULL;
  Lstrict = false;
  Lkmc = false;
  Ladapt = false;
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
    } else if (strcmp(arg[iarg],"adapt") == 0) {
      if (iarg+2 > narg) error->all("Illegal sweep_style command");
      if (strcmp(arg[iarg+1],"yes") == 0) Ladapt = true;
      else if (strcmp(arg[iarg+1],"no") == 0) Ladapt = false;
      else error->all("Illegal sweep_style command");
      iarg += 2;
    } else error->all("Illegal sweep_style command");
  }

  // set sweep method function ptr

  if (Lkmc && (Lmask || Lstrict))
      error->all("Invalid combination of sweep flags");
  if (Ladapt && !Lkmc)
      error->all("Invalid combination of sweep flags");

  if (Lkmc) sweeper = &SweepLattice::sweep_sector_kmc;
  else {
    if (!Lmask && !Lstrict) sweeper = &SweepLattice::sweep_sector;
    if (Lmask && !Lstrict) sweeper = &SweepLattice::sweep_sector_mask;
    if (!Lmask && Lstrict) sweeper = &SweepLattice::sweep_sector_strict;
    if (Lmask && Lstrict) sweeper = &SweepLattice::sweep_sector_mask_strict;
  }

  // initializations

  int nsectormax = 8;
  for (int i = 0; i < nsectormax; i++) {
    sector[i].solve = NULL;
    sector[i].site2i = NULL;
    sector[i].i2site = NULL;
    sector[i].propensity = NULL;
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
    memory->sfree(sector[isector].site2i);
    memory->sfree(sector[isector].i2site);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->sfree(sector[isector].border);
  }
}

/* ---------------------------------------------------------------------- */

void SweepLattice::init()
{
  int i,j,k,m;

  applattice = (AppLattice *) app;

  // error check

  if (Lmask && applattice->temperature > 0.0)
    error->all("Cannot mask sweeping with non-zero temperature");

  // get application params

  dimension = applattice->dimension;
  nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  int *id = applattice->id;
  double **xyz = applattice->xyz;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  int delpropensity = applattice->delpropensity;
  int delevent = applattice->delevent;

  delcol = delpropensity+1;
  if (Lstrict) ncolor = delcol*delcol;
  else ncolor = 1;

  // setup sectors

  if (dimension == 2) nsector = 4;
  else if (dimension == 3) nsector = 8;

  for (int isector = 0; isector < nsector; isector++) {
    sector[isector].nlocal = 0;
    sector[isector].nmax = 0;
    delete sector[isector].solve;
    memory->sfree(sector[isector].site2i);
    memory->sfree(sector[isector].i2site);
    memory->sfree(sector[isector].propensity);
    memory->sfree(sector[isector].sites);
    memory->sfree(sector[isector].border);
    sector[isector].solve = NULL;
    sector[isector].propensity = NULL;
    sector[isector].site2i = NULL;
    sector[isector].i2site = NULL;
    sector[isector].sites = NULL;
    sector[isector].border = NULL;
  }

  int iwhich,jwhich,kwhich,isector;
  double xmid = 0.5 * (applattice->subxlo + applattice->subxhi);
  double ymid = 0.5 * (applattice->subylo + applattice->subyhi);
  double zmid = 0.5 * (applattice->subzlo + applattice->subzhi);

  // assign each owned site to a sector (4/8 in 2d/3d)
  // setup site2i for each sector

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
				 "sweep:site2i");
    }
    sector[isector].site2i[sector[isector].nlocal++] = i;
  }

  // setup i2site for each sector
  // 0 to nsite-1 for owned points in sector, else -1

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++) {
      int nsites = sector[isector].nlocal;

      sector[isector].i2site = 
	(int *) memory->smalloc(ntotal*sizeof(int),"sweep:i2site");
      for (i = 0; i < ntotal; i++) sector[isector].i2site[i] = -1;
      for (m = 0; m < nsites; m++) 
	sector[isector].i2site[sector[isector].site2i[m]] = m;
    }
  }

  // setup mask array for owned and ghost values
  // ghost values written to in app::site_event_rejection

  if (Lmask && mask == NULL) {
    mask = (char *) memory->smalloc(ntotal*sizeof(char),"sweep:mask");
    for (i = 0; i < nlocal; i++) mask[i] = 0;
  }

  // setup border array used by masking and KMC solver
  // nborder = # of sites in sector influenced by site outside sector
  // border = list of these sites, stored as lattice indices

  if (Lmask || Lkmc) {
    for (int isector = 0; isector < nsector; isector++)
      sector[isector].nborder = find_border_sites(isector,numneigh,neighbor);
  }

  // allocate sites array used by KMC solver for border sites

  if (Lkmc) {
    for (int isector = 0; isector < nsector; isector++)
      sector[isector].sites =
	(int *) memory->smalloc(sector[isector].nborder*sizeof(int),
				"sweep:sites");
  }

  // setup one RNG per site
  // only owned values referenced in strict() methods

  if (Lstrict && ranlat == NULL) {
    ranlat = (RandomPark *) memory->smalloc(nlocal*sizeof(RandomPark *),
					    "sweep:ranlat");
    for (i = 0; i < nlocal; i++) ranlat[i].init(seed+id[i]);
  }

  // init communication for ghost sites, requires site2i

  comm->init(this,delpropensity,delevent,NULL);

  // setup KMC solver and propensity arrays for each sector
  // propensity init requires all ghost cell info via comm->all()

  if (Lkmc) {
    comm->all();

    for (int isector = 0; isector < nsector; isector++) {
      sector[isector].solve = solve->clone();

      int nsites = sector[isector].nlocal;
      sector[isector].propensity = 
	(double *) memory->smalloc(nsites*sizeof(double),"sweep:propensity");
      for (m = 0; m < nsites; m++)
	sector[isector].propensity[m] =
	  applattice->site_propensity(sector[isector].site2i[m]);

      sector[isector].solve->init(nsites,sector[isector].propensity);
    }
  }

  // compute deln0, which controls future timesteps

  if (Lkmc && Ladapt) {
    pmax = 0.0;
    int ntmp ;
    double ptmp;
    for (int isector = 0; isector < nsector; isector++) {
      ntmp = sector[isector].solve->get_num_active();
      if (ntmp > 0) {
	ptmp = sector[isector].solve->get_total_propensity();
	ptmp /= ntmp;
	pmax = MAX(ptmp,pmax);
      }
    }
    MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    deln0 = pmaxall*delt;
  }
}

/* ----------------------------------------------------------------------
   perform single sweep over entire domain by sectors
   interleave communication
 ------------------------------------------------------------------------- */

void SweepLattice::do_sweep(double &dt)
{
  if (Ladapt) pmax = 0.0;
  if (!Lkmc) applattice->ntimestep++;

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

  // adjust KMC threshold time

  if (Ladapt) {
    MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_SUM,world);
    if (pmaxall > 0.0) {
      delt = deln0/pmaxall;
    }
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector(int icolor, int isector)
{
  int i,m;

  int *site2i = sector[isector].site2i;
  int nlocal = sector[isector].nlocal;

  for (m = 0; m < nlocal; m++) {
    i = site2i[m];
    applattice->site_event_rejection(i,random);
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites with masking
   performs events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_mask(int icolor, int isector)
{
  int i,m;
  
  boundary_clear_mask(isector);
  
  int *site2i = sector[isector].site2i;
  int nlocal = sector[isector].nlocal;

  for (m = 0; m < nlocal; m++) {
    i = site2i[m];
    if (mask[i]) continue;
    applattice->site_event_rejection(i,random);
  }
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color
   perform events with rejection
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_strict(int icolor, int isector)
{
  error->all("Sweep option not yet supported");
}

/* ----------------------------------------------------------------------
   sweep over one sector of sites of one color with masking
   perform events with rejection
   update mask values
 ------------------------------------------------------------------------- */
   
void SweepLattice::sweep_sector_mask_strict(int icolor, int isector)
{
  error->all("Sweep option not yet supported");
}

/* ----------------------------------------------------------------------
   generate events in each sector using KMC solver
 ------------------------------------------------------------------------- */

void SweepLattice::sweep_sector_kmc(int icolor, int isector)
{
  double dt,time;
  int done,isite,i;

  // extract sector specific info from struct

  int *border = sector[isector].border;
  int nborder = sector[isector].nborder;

  Solve *solve = sector[isector].solve;
  double *propensity = sector[isector].propensity;
  int *site2i = sector[isector].site2i;
  int *i2site = sector[isector].i2site;
  int *sites = sector[isector].sites;

  // temporarily reset values in applattice
  // solver, propensity, i2site array

  Solve *hold_solve = applattice->solve;
  applattice->solve = solve;
  applattice->propensity = propensity;
  applattice->i2site = i2site;

  // update propensities for sites which neighbor a site outside sector
  // necessary since outside sites may have changed
  // attribute this chunk of time to comm, b/c due to decomposition

  timer->stamp();

  int nsites = 0;

  for (int m = 0; m < nborder; m++) {
    i = border[m];
    isite = i2site[i];
    sites[nsites++] = isite;
    propensity[isite] = applattice->site_propensity(i);
  }

  solve->update(nsites,sites,propensity);
  timer->stamp(TIME_COMM);

  // execute events until time threshhold reached

  done = 0;
  time = 0.0;
  while (!done) {
    timer->stamp();
    isite = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    if (isite < 0) done = 1;
    else {
      time += dt;	
      if (time >= delt) {
    	done = 1;
      } else {
	i = site2i[isite];
	applattice->site_event(i,random);
	applattice->ntimestep++;
      }
      timer->stamp(TIME_APP);
    }
  }
 
  // compute deln0, which controls future timesteps

  if (Ladapt) {
    int ntmp = applattice->solve->get_num_active();
    if (ntmp > 0) {
      double ptmp = applattice->solve->get_total_propensity();
      ptmp /= ntmp;
      pmax = MAX(ptmp,pmax);
    }
  }

  // restore applattice values

  applattice->solve = hold_solve;
  applattice->propensity = NULL;
  applattice->i2site = NULL;
}

/* ----------------------------------------------------------------------
   create list of border sites for a sector
   border site = site in sector with a neighbor outside the sector
   neighbor can be another owned site (outside sector) or a ghost
   border = lattice index of the sites
 ------------------------------------------------------------------------- */

int SweepLattice::find_border_sites(int isector,
				    int *numneigh, int **neighbor)
{
  int nsites = sector[isector].nlocal;
  int *site2i = sector[isector].site2i;

  int *border = NULL;
  int nborder = 0;
  int nmax = 0;

  // put lattice index of each site into hash

  int i,j,m;

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;
  for (int m = 0; m < nsites; m++)
    hash.insert(std::pair<int,int> (site2i[m],0));
  
  // add site to border if has neighbor not in hash

  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    for (j = 0; j < numneigh[i]; j++)
      if (hash.find(neighbor[i][j]) == hash.end()) break;
    if (j == numneigh[i]) continue;
    
    if (nborder == nmax) {
      nmax += DELTA;
      border = (int *) memory->srealloc(border,nmax*sizeof(int),
					"sweep:border");
    }
    border[nborder++] = i;
  }

  sector[isector].border = border;
  return nborder;
}
  
/* ----------------------------------------------------------------------
   unset all mask values of owned sites in isector whose propensity
     could change due to events on sites one neighbor outside the sector
   border list stores indices of these sites
   their mask value may be out of date, due to state change in other sectors
 ------------------------------------------------------------------------- */

void SweepLattice::boundary_clear_mask(int isector)
{
  int *border = sector[isector].border;
  int nborder = sector[isector].nborder;

  for (int m = 0; m < nborder; m++) mask[border[m]] = 0;
}
