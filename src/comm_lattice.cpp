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

#include "spktype.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice.h"
#include "app.h"
#include "app_lattice.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

#define DELTA 16384

/* ---------------------------------------------------------------------- */

CommLattice::CommLattice(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  allswap = NULL;
  reverseswap = NULL;

  nsector = 0;
  sectorswap = NULL;
  sectorreverseswap = NULL;
}

/* ---------------------------------------------------------------------- */

CommLattice::~CommLattice()
{
  if (allswap) free_swap(allswap);
  if (reverseswap) free_swap(reverseswap);
  if (sectorswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
    delete [] sectorswap;
  }
  if (sectorreverseswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorreverseswap[i]);
    delete [] sectorreverseswap;
  }
}

/* ----------------------------------------------------------------------
   setup comm pattern
   3 kinds of patterns: all ghosts, sector ghosts, sector reverse
   if nsector_request = 1, just ghosts for entire proc domain 
   if nsector_request > 1, do all and sector ghosts, reverse if needed
   array = NULL = communicate iarray/darray from app
   array = non-NULL = communicate passed-in array (from diagnostic)
------------------------------------------------------------------------- */

void CommLattice::init(int nsector_request, int delpropensity, int delevent,
		       int *array) 
{
  delghost = delpropensity;
  delreverse = delevent;

  AppLattice *applattice = (AppLattice *) app;

  ninteger = app->ninteger;
  ndouble = app->ndouble;
  iarray = app->iarray;
  darray = app->darray;

  if (array) {
    site_only = 1;
    site = array;
  } else if (ninteger == 1 && ndouble == 0) {
    site_only = 1;
    site = iarray[0];
  } else site_only = 0;

  // clear out old swaps

  if (allswap) free_swap(allswap);
  if (reverseswap) free_swap(reverseswap);
  if (sectorswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
    delete [] sectorswap;
  }
  if (sectorreverseswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorreverseswap[i]);
    delete [] sectorreverseswap;
  }
  allswap = NULL;
  reverseswap = NULL;
  sectorswap = NULL;
  sectorreverseswap = NULL;

  // create new swaps as requested

  allswap = create_swap_all();
  reverseswap = create_swap_all_reverse();

  nsector = nsector_request;
  if (nsector > 1) {
    sectorswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorswap[i] = create_swap_sector(applattice->set[i].nlocal,
					 applattice->set[i].site2i);
  }
  if (delreverse && nsector > 1) {
    sectorreverseswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorreverseswap[i] = 
	create_swap_sector_reverse(applattice->set[i].nlocal,
				   applattice->set[i].site2i);
  }
}

/* ----------------------------------------------------------------------
   acquire ghost values for entire proc sub-domain
------------------------------------------------------------------------- */

void CommLattice::all()
{
  if (site_only) perform_swap_site(allswap);
  else if (ndouble == 0) perform_swap_int(allswap);
  else if (ninteger == 0) perform_swap_double(allswap);
  else perform_swap_general(allswap);
}

/* ----------------------------------------------------------------------
   reverse communicate changed border values for entire proc sub-domain
------------------------------------------------------------------------- */

void CommLattice::all_reverse()
{
  if (delreverse == 0) return;

  if (site_only) perform_swap_site(reverseswap);
  else if (ndouble == 0) perform_swap_int(reverseswap);
  else if (ninteger == 0) perform_swap_double(reverseswap);
  else perform_swap_general(reverseswap);
}

/* ----------------------------------------------------------------------
   acquire ghost values for one sector
------------------------------------------------------------------------- */

void CommLattice::sector(int isector)
{
  if (site_only) perform_swap_site(sectorswap[isector]);
  else if (ndouble == 0) perform_swap_int(sectorswap[isector]);
  else if (ninteger == 0) perform_swap_double(sectorswap[isector]);
  else perform_swap_general(sectorswap[isector]);
}

/* ----------------------------------------------------------------------
   reverse communicate changed border values for one sector
------------------------------------------------------------------------- */

void CommLattice::reverse_sector(int isector)
{
  if (delreverse == 0) return;

  if (site_only) perform_swap_site(sectorreverseswap[isector]);
  else if (ndouble == 0) perform_swap_int(sectorreverseswap[isector]);
  else if (ninteger == 0) perform_swap_double(sectorreverseswap[isector]);
  else perform_swap_general(sectorreverseswap[isector]);
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost sites for my entire subdomain
------------------------------------------------------------------------- */

CommLattice::Swap *CommLattice::create_swap_all()
{
  int i;

  AppLattice *applattice = (AppLattice *) app;
  int nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  tagint *id = app->id;

  // buf = list of global IDs I need to receive

  Site *buf = (Site *) memory->smalloc(nghost*sizeof(Site),"comm:buf");
  int nsite = 0;

  for (i = nlocal; i < ntotal; i++) {
    buf[nsite].id_global = id[i];
    buf[nsite].index_local = i;
    buf[nsite].proc = -1;
    nsite++;
  }

  // grow buf to size of max ghosts on any proc

  int maxsite;
  MPI_Allreduce(&nsite,&maxsite,1,MPI_INT,MPI_MAX,world);

  buf = (Site *) memory->srealloc(buf,maxsite*sizeof(Site),"comm:buf");

  // create swap based on list of recvs

  Swap *swap = new Swap;

  create_send_from_recv(nsite,maxsite,buf,swap);
  create_recv_from_list(nsite,buf,swap);

  memory->sfree(buf);

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   reverse comm for my entire subdomain
------------------------------------------------------------------------- */

CommLattice::Swap *CommLattice::create_swap_all_reverse()
{
  int i,j;

  AppLattice *applattice = (AppLattice *) app;
  int nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  tagint *id = app->id;
  int *owner = applattice->owner;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  // flag ghost sites with -1
  // flag owned sites with 0

  int *flag;
  memory->create(flag,ntotal,"comm:flag");
  for (i = 0; i < ntotal; i++) flag[i] = -1;
  for (i = 0; i < nlocal; i++) flag[i] = 0;

  // flag ghost sites up to delreverse hops from subdomain with positive number
  // assumes ghost sites have a neighbor list

  for (int ilayer = 0; ilayer < delreverse; ilayer++) {
    for (i = 0; i < ntotal; i++) {
      if (flag[i] != ilayer) continue;
      for (j = 0; j < numneigh[i]; j++) {
	if (neighbor[i][j] < nlocal) continue;
	if (flag[neighbor[i][j]] == -1) flag[neighbor[i][j]] = ilayer+1;
      }
    }
  }

  // buf = list of global IDs I will send

  Site *buf = (Site *) memory->smalloc(nghost*sizeof(Site),"comm:buf");
  int nsite = 0;

  for (i = nlocal; i < ntotal; i++) {
    if (flag[i] <= 0) continue;
    buf[nsite].id_global = id[i];
    buf[nsite].index_local = i;
    buf[nsite].proc = owner[i];
    nsite++;
  }

  memory->destroy(flag);

  // grow buf to size of max sites on any proc

  int maxsite;
  MPI_Allreduce(&nsite,&maxsite,1,MPI_INT,MPI_MAX,world);

  buf = (Site *) memory->srealloc(buf,maxsite*sizeof(Site),"comm:buf");

  // create swap based on list of sends

  Swap *swap = new Swap;

  create_send_from_list(nsite,buf,swap);
  create_recv_from_send(nsite,maxsite,buf,swap);

  memory->sfree(buf);

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost sites for a single sector of my subdomain
------------------------------------------------------------------------- */

CommLattice::Swap *CommLattice::create_swap_sector(int nsites, int *site2i)
{
  int i,j,m;

  AppLattice *applattice = (AppLattice *) app;
  int nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  tagint *id = app->id;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  // flag sites with -1 that are not in sector
  // flag sites with 0 that are in sector

  int *flag;
  memory->create(flag,ntotal,"comm:flag");
  for (i = 0; i < ntotal; i++) flag[i] = -1;
  for (m = 0; m < nsites; m++) flag[site2i[m]] = 0;

  // flag ghost sites up to delghost hops from sector with positive number
  // assumes ghost sites have a neighbor list

  for (int ilayer = 0; ilayer < delghost; ilayer++) {
    for (i = 0; i < ntotal; i++) {
      if (flag[i] != ilayer) continue;
      for (j = 0; j < numneigh[i]; j++) {
	if (neighbor[i][j] < nlocal) continue;
	if (flag[neighbor[i][j]] == -1) flag[neighbor[i][j]] = ilayer+1;
      }
    }
  }

  // buf = list of global IDs I need to receive

  Site *buf = (Site *) memory->smalloc(nghost*sizeof(Site),"comm:buf");
  int nsite = 0;

  for (i = nlocal; i < ntotal; i++) {
    if (flag[i] <= 0) continue;
    buf[nsite].id_global = id[i];
    buf[nsite].index_local = i;
    buf[nsite].proc = -1;
    nsite++;
  }

  memory->destroy(flag);

  // grow buf to size of max sites on any proc

  int maxsite;
  MPI_Allreduce(&nsite,&maxsite,1,MPI_INT,MPI_MAX,world);

  buf = (Site *) memory->srealloc(buf,maxsite*sizeof(Site),"comm:buf");

  // create swap based on list of recvs

  Swap *swap = new Swap;

  create_send_from_recv(nsite,maxsite,buf,swap);
  create_recv_from_list(nsite,buf,swap);

  memory->sfree(buf);

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   reverse comm of a single sector
------------------------------------------------------------------------- */

CommLattice::Swap *CommLattice::create_swap_sector_reverse(int nsites,
							   int *site2i)
{
  int i,j,m;

  AppLattice *applattice = (AppLattice *) app;
  int nlocal = applattice->nlocal;
  int nghost = applattice->nghost;
  int ntotal = nlocal + nghost;

  tagint *id = app->id;
  int *owner = applattice->owner;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  // flag sites with -1 that are not in sector
  // flag sites with 0 that are in sector

  int *flag;
  memory->create(flag,ntotal,"comm:flag");
  for (i = 0; i < ntotal; i++) flag[i] = -1;
  for (m = 0; m < nsites; m++) flag[site2i[m]] = 0;

  // flag ghost sites up to delreverse hops from sector with positive number
  // assumes ghost sites have a neighbor list

  for (int ilayer = 0; ilayer < delreverse; ilayer++) {
    for (i = 0; i < ntotal; i++) {
      if (flag[i] != ilayer) continue;
      for (j = 0; j < numneigh[i]; j++) {
	if (neighbor[i][j] < nlocal) continue;
	if (flag[neighbor[i][j]] == -1) flag[neighbor[i][j]] = ilayer+1;
      }
    }
  }

  // buf = list of global IDs I will send

  Site *buf = (Site *) memory->smalloc(nghost*sizeof(Site),"comm:buf");
  int nsite = 0;

  for (i = nlocal; i < ntotal; i++) {
    if (flag[i] <= 0) continue;
    buf[nsite].id_global = id[i];
    buf[nsite].index_local = i;
    buf[nsite].proc = owner[i];
    nsite++;
  }

  memory->destroy(flag);

  // grow buf to size of max sites on any proc

  int maxsite;
  MPI_Allreduce(&nsite,&maxsite,1,MPI_INT,MPI_MAX,world);

  buf = (Site *) memory->srealloc(buf,maxsite*sizeof(Site),"comm:buf");

  // create swap based on list of sends

  Swap *swap = new Swap;

  create_send_from_list(nsite,buf,swap);
  create_recv_from_send(nsite,maxsite,buf,swap);

  memory->sfree(buf);

  return swap;
}

/* ----------------------------------------------------------------------
   create send portion of a Swap communication pattern
   start with list of sites I need to recv
   circulate list to all procs
   when comes back to me, will also be flagged with who I recv sites from
------------------------------------------------------------------------- */

void CommLattice::create_send_from_recv(int nsite, int maxsite,
					Site *buf, Swap *swap)
{
  int i,isend;

  AppLattice *applattice = (AppLattice *) app;
  tagint *id = app->id;
  int nlocal = applattice->nlocal;

  Site *bufcopy = (Site *) 
    memory->smalloc(maxsite*sizeof(Site),"comm:bufcopy");

  // put my entire list of owned site IDs in a hash

  std::map<tagint,int>::iterator loc;
  std::map<tagint,int> hash;

  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<tagint,int> (id[i],i));

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // allocate send lists

  int nsend = 0;
  int *sproc = new int[nprocs];
  int *scount = new int[nprocs];
  int *smax = new int[nprocs];
  int **sindex = new int*[nprocs];
  for (int i = 0; i < nprocs; i++) {
    smax[i] = 0;
    sindex[i] = NULL;
  }

  // cycle request list around ring of procs back to self
  // when receive it, add sites I own to my send lists, also fill in proc field
  // original_proc = proc to send info to during swap

  MPI_Request request;
  MPI_Status status;

  int size = nsite;
  int original_proc = me;

  for (int loop = 0; loop < nprocs; loop++) {
    original_proc--;
    if (original_proc < 0) original_proc = nprocs-1;

    if (me != next) {
      MPI_Irecv(bufcopy,maxsite*sizeof(Site),MPI_CHAR,prev,0,world,&request);
      MPI_Send(buf,size*sizeof(Site),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Site);
      memcpy(buf,bufcopy,size*sizeof(Site));
    }

    for (i = 0; i < size; i++) {
      if (buf[i].proc >= 0) continue;
      loc = hash.find(buf[i].id_global);
      if (loc == hash.end()) continue;
      buf[i].proc = me;

      // if original_proc not last one in sproc[], add proc to send list

      if (nsend == 0 || sproc[nsend-1] != original_proc) {
	sproc[nsend] = original_proc;
	scount[nsend] = 0;
	sindex[nsend] = NULL;
	nsend++;
      }

      // add index to send list going to a particular proc
      // grow sindex array if necessary

      isend = nsend - 1;
      if (scount[isend] == smax[isend]) {
	smax[isend] += DELTA;
	memory->grow(sindex[isend],smax[isend],"comm:sindex");
      }
      sindex[isend][scount[isend]] = loc->second;
      scount[isend]++;
    }
  }

  // allocate sbuf
  // sizes depend on number of ints and doubles stored per site

  int max = 0;
  for (i = 0; i < nsend; i++) max = MAX(max,scount[i]);

  int *sibuf = NULL;
  double *sdbuf = NULL;
  if (max) {
    if (site_only) sibuf = new int[max];
    else if (ndouble == 0) sibuf = new int[ninteger*max];
    else if (ninteger == 0) sdbuf = new double[ndouble*max];
    else sdbuf = new double[(ninteger+ndouble)*max];
  }

  // clean up

  memory->sfree(bufcopy);

  // fill in swap data structure

  swap->nsend = nsend;
  swap->sproc = sproc;
  swap->scount = scount;
  swap->smax = smax;
  swap->sindex = sindex;
  swap->sibuf = sibuf;
  swap->sdbuf = sdbuf;
}

/* ----------------------------------------------------------------------
   create send portion of a Swap communication pattern
   create from list of sites I send and procs who own them
------------------------------------------------------------------------- */

void CommLattice::create_send_from_list(int nsite, Site *buf, Swap *swap)
{
  int i,isend;
  std::map<int,int>::iterator loc;
  std::map<int,int> hash;

  // allocate send lists

  int nsend = 0;
  int *sproc = new int[nprocs];
  int *scount = new int[nprocs];
  int *smax = new int[nprocs];
  int **sindex = new int*[nprocs];
  for (int i = 0; i < nprocs; i++) {
    smax[i] = 0;
    sindex[i] = NULL;
  }

  // extract info for send lists

  for (i = 0; i < nsite; i++) {

    // if sending proc not in sproc[], add proc to send list and to hash
    // isend = location of this proc in send lists

    loc = hash.find(buf[i].proc);
    if (loc == hash.end()) {
      sproc[nsend] = buf[i].proc;
      scount[nsend] = 0;
      sindex[nsend] = NULL;
      hash.insert(std::pair<int,int> (buf[i].proc,nsend));
      isend = nsend;
      nsend++;
    } else isend = loc->second;
    
    // add index to send list going to a particular proc
    // grow rindex array if necessary

    if (scount[isend] == smax[isend]) {
      smax[isend] += DELTA;
      memory->grow(sindex[isend],smax[isend],"comm:sindex");
    }
    sindex[isend][scount[isend]] = buf[i].index_local;
    scount[isend]++;
  }

  // allocate sbuf
  // sizes depend on number of ints and doubles stored per site

  int max = 0;
  for (i = 0; i < nsend; i++) max = MAX(max,scount[i]);

  int *sibuf = NULL;
  double *sdbuf = NULL;
  if (max) {
    if (site_only) sibuf = new int[max];
    else if (ndouble == 0) sibuf = new int[ninteger*max];
    else if (ninteger == 0) sdbuf = new double[ndouble*max];
    else sdbuf = new double[(ninteger+ndouble)*max];
  }

  // fill in swap data structure

  swap->nsend = nsend;
  swap->sproc = sproc;
  swap->scount = scount;
  swap->smax = smax;
  swap->sindex = sindex;
  swap->sibuf = sibuf;
  swap->sdbuf = sdbuf;
}

/* ----------------------------------------------------------------------
   create recv portion of a Swap communication pattern
   start with list of sites I will send
   circulate list to all procs
------------------------------------------------------------------------- */

void CommLattice::create_recv_from_send(int nsite, int maxsite,
					Site *buf, Swap *swap)
{
  int i,irecv;

  AppLattice *applattice = (AppLattice *) app;
  tagint *id = app->id;
  int nlocal = applattice->nlocal;

  Site *bufcopy = (Site *) 
    memory->smalloc(maxsite*sizeof(Site),"comm:bufcopy");

  // put my entire list of owned site IDs in a hash

  std::map<tagint,int>::iterator loc;
  std::map<tagint,int> hash;

  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<tagint,int> (id[i],i));

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // allocate recv lists

  int nrecv = 0;
  int *rproc = new int[nprocs];
  int *rcount = new int[nprocs];
  int *rmax = new int[nprocs];
  int **rindex = new int*[nprocs];
  int **ribuf = new int*[nprocs];
  double **rdbuf = new double*[nprocs];
  for (int i = 0; i < nprocs; i++) {
    rmax[i] = 0;
    rindex[i] = NULL;
    ribuf[i] = NULL;
    rdbuf[i] = NULL;
  }

  // cycle send list around ring of procs back to self
  // when receive it, add sites I own to my recv lists
  // original_proc = proc to recv info from during swap

  MPI_Request request;
  MPI_Status status;

  int size = nsite;
  int original_proc = me;

  for (int loop = 0; loop < nprocs; loop++) {
    original_proc--;
    if (original_proc < 0) original_proc = nprocs-1;

    if (me != next) {
      MPI_Irecv(bufcopy,maxsite*sizeof(Site),MPI_CHAR,prev,0,world,&request);
      MPI_Send(buf,size*sizeof(Site),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Site);
      memcpy(buf,bufcopy,size*sizeof(Site));
    }

    for (i = 0; i < size; i++) {
      loc = hash.find(buf[i].id_global);
      if (loc == hash.end()) continue;

      // if original_proc not last one in rproc[], add proc to recv list

      if (nrecv == 0 || rproc[nrecv-1] != original_proc) {
	rproc[nrecv] = original_proc;
	rcount[nrecv] = 0;
	rindex[nrecv] = NULL;
	nrecv++;
      }

      // add index to recv list from a particular proc
      // grow rindex array if necessary
      
      irecv = nrecv - 1;
      if (rcount[irecv] == rmax[irecv]) {
	rmax[irecv] += DELTA;
	memory->grow(rindex[irecv],rmax[irecv],"comm:rindex");
      }
      rindex[irecv][rcount[irecv]] = loc->second;
      rcount[irecv]++;
    }
  }

  // allocate rbuf[i]
  // sizes depend on number of ints and doubles stored per site

  for (i = 0; i < nrecv; i++) {
    if (site_only) ribuf[i] = new int[rcount[i]];
    else if (ndouble == 0) ribuf[i] = new int[ninteger*rcount[i]];
    else if (ninteger == 0) rdbuf[i] = new double[ndouble*rcount[i]];
    else rdbuf[i] = new double[(ninteger+ndouble)*rcount[i]];
  }

  // clean up

  memory->sfree(bufcopy);

  // fill in swap data struct

  swap->nrecv = nrecv;
  swap->rproc = rproc;
  swap->rcount = rcount;
  swap->rmax = rmax;
  swap->rindex = rindex;
  swap->ribuf = ribuf;
  swap->rdbuf = rdbuf;

  if (swap->nrecv) {
    swap->request = new MPI_Request[swap->nrecv];
    swap->status = new MPI_Status[swap->nrecv];
  } else {
    swap->request = NULL;
    swap->status = NULL;
  }
}

/* ----------------------------------------------------------------------
   create recv portion of a Swap communication pattern
   create from list of sites I need to recv and procs who will send them
   no communication required
------------------------------------------------------------------------- */

void CommLattice::create_recv_from_list(int nsite, Site *buf, Swap *swap)
{
  int i,irecv;
  std::map<int,int>::iterator loc;
  std::map<int,int> hash;

  // allocate recv lists

  int nrecv = 0;
  int *rproc = new int[nprocs];
  int *rcount = new int[nprocs];
  int *rmax = new int[nprocs];
  int **rindex = new int*[nprocs];
  int **ribuf = new int*[nprocs];
  double **rdbuf = new double*[nprocs];
  for (int i = 0; i < nprocs; i++) {
    rmax[i] = 0;
    rindex[i] = NULL;
    ribuf[i] = NULL;
    rdbuf[i] = NULL;
  }

  // extract info for recv lists
  // error if any site is not filled in

  for (i = 0; i < nsite; i++) {
    if (buf[i].proc == -1) 
      error->one(FLERR,"Site-site interaction was not found");

    // if sending proc not in rproc[], add proc to recv list and to hash
    // irecv = location of this proc in recv lists

    loc = hash.find(buf[i].proc);
    if (loc == hash.end()) {
      rproc[nrecv] = buf[i].proc;
      rcount[nrecv] = 0;
      rindex[nrecv] = NULL;
      hash.insert(std::pair<int,int> (buf[i].proc,nrecv));
      irecv = nrecv;
      nrecv++;
    } else irecv = loc->second;
    
    // add index to recv list from a particular proc
    // grow rindex array if necessary

    if (rcount[irecv] == rmax[irecv]) {
      rmax[irecv] += DELTA;
      memory->grow(rindex[irecv],rmax[irecv],"comm:rindex");
    }
    rindex[irecv][rcount[irecv]] = buf[i].index_local;
    rcount[irecv]++;
  }

  // allocate rbuf[i]
  // sizes depend on number of ints and doubles stored per site

  for (i = 0; i < nrecv; i++) {
    if (site_only) ribuf[i] = new int[rcount[i]];
    else if (ndouble == 0) ribuf[i] = new int[ninteger*rcount[i]];
    else if (ninteger == 0) rdbuf[i] = new double[ndouble*rcount[i]];
    else rdbuf[i] = new double[(ninteger+ndouble)*rcount[i]];
  }

  // fill in swap data struct

  swap->nrecv = nrecv;
  swap->rproc = rproc;
  swap->rcount = rcount;
  swap->rmax = rmax;
  swap->rindex = rindex;
  swap->ribuf = ribuf;
  swap->rdbuf = rdbuf;

  if (swap->nrecv) {
    swap->request = new MPI_Request[swap->nrecv];
    swap->status = new MPI_Status[swap->nrecv];
  } else {
    swap->request = NULL;
    swap->status = NULL;
  }
}

/* ----------------------------------------------------------------------
   free a Swap object
------------------------------------------------------------------------- */

void CommLattice::free_swap(Swap *swap)
{
  delete [] swap->sproc;
  delete [] swap->scount;
  delete [] swap->smax;
  for (int i = 0; i < swap->nsend; i++) memory->destroy(swap->sindex[i]);
  delete [] swap->sindex;
  delete [] swap->sibuf;
  delete [] swap->sdbuf;

  delete [] swap->rproc;
  delete [] swap->rcount;
  delete [] swap->rmax;
  for (int i = 0; i < swap->nrecv; i++) memory->destroy(swap->rindex[i]);
  delete [] swap->rindex;
  for (int i = 0; i < swap->nrecv; i++) delete [] swap->ribuf[i];
  for (int i = 0; i < swap->nrecv; i++) delete [] swap->rdbuf[i];
  delete [] swap->ribuf;
  delete [] swap->rdbuf;

  delete [] swap->request;
  delete [] swap->status;

  delete swap;
}

/* ----------------------------------------------------------------------
   communicate site values via Swap instructions
   use site array = iarray[0] as source/destination
------------------------------------------------------------------------- */

void CommLattice::perform_swap_site(Swap *swap)
{
  int i,j;
  int *index;
  int *buf;

  // post receives

  for (i = 0; i < swap->nrecv; i++)
    MPI_Irecv(swap->ribuf[i],swap->rcount[i],MPI_INT,swap->rproc[i],0,world,
	      &swap->request[i]);

  // pack data to send to each proc and send it

  for (i = 0; i < swap->nsend; i++) {
    index = swap->sindex[i];
    buf = swap->sibuf;
    for (j = 0; j < swap->scount[i]; j++)
      buf[j] = site[index[j]];
    MPI_Send(buf,swap->scount[i],MPI_INT,swap->sproc[i],0,world);
  }

  // wait on incoming messages

  if (swap->nrecv) MPI_Waitall(swap->nrecv,swap->request,swap->status);

  // unpack received buffers of data from each proc

  for (i = 0; i < swap->nrecv; i++) {
    index = swap->rindex[i];
    buf = swap->ribuf[i];
    for (j = 0; j < swap->rcount[i]; j++)
      site[index[j]] = buf[j];
  }
}

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use iarray as source/destination, all integer data
------------------------------------------------------------------------- */

void CommLattice::perform_swap_int(Swap *swap)
{
  int i,j,m,n;
  int *index,*buf,*vector;

  // post receives

  for (i = 0; i < swap->nrecv; i++)
    MPI_Irecv(swap->ribuf[i],ninteger*swap->rcount[i],MPI_INT,
	      swap->rproc[i],0,world,&swap->request[i]);

  // pack data to send to each proc and send it

  for (i = 0; i < swap->nsend; i++) {
    index = swap->sindex[i];
    buf = swap->sibuf;
    m = 0;
    for (n = 0; n < ninteger; n++) {
      vector = iarray[n];
      for (j = 0; j < swap->scount[i]; j++) buf[m++] = vector[index[j]];
    }
    MPI_Send(buf,ninteger*swap->scount[i],MPI_INT,swap->sproc[i],0,world);
  }

  // wait on incoming messages

  if (swap->nrecv) MPI_Waitall(swap->nrecv,swap->request,swap->status);

  // unpack received buffers of data from each proc

  for (i = 0; i < swap->nrecv; i++) {
    index = swap->rindex[i];
    buf = swap->ribuf[i];
    m = 0;
    for (n = 0; n < ninteger; n++) {
      vector = iarray[n];
      for (j = 0; j < swap->rcount[i]; j++) vector[index[j]] = buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use darray as source/destination, all double data
------------------------------------------------------------------------- */

void CommLattice::perform_swap_double(Swap *swap)
{
  int i,j,m,n;
  int *index;
  double *buf,*vector;

  // post receives

  for (i = 0; i < swap->nrecv; i++)
    MPI_Irecv(swap->rdbuf[i],ndouble*swap->rcount[i],MPI_DOUBLE,
	      swap->rproc[i],0,world,&swap->request[i]);

  // pack data to send to each proc and send it

  for (i = 0; i < swap->nsend; i++) {
    index = swap->sindex[i];
    buf = swap->sdbuf;
    m = 0;
    for (n = 0; n < ndouble; n++) {
      vector = darray[n];
      for (j = 0; j < swap->scount[i]; j++) buf[m++] = vector[index[j]];
    }
    MPI_Send(buf,ndouble*swap->scount[i],MPI_DOUBLE,swap->sproc[i],0,world);
  }

  // wait on incoming messages

  if (swap->nrecv) MPI_Waitall(swap->nrecv,swap->request,swap->status);

  // unpack received buffers of data from each proc

  for (i = 0; i < swap->nrecv; i++) {
    index = swap->rindex[i];
    buf = swap->rdbuf[i];
    m = 0;
    for (n = 0; n < ndouble; n++) {
      vector = darray[n];
      for (j = 0; j < swap->rcount[i]; j++) vector[index[j]] = buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
   use iarray and darray as source/destination, mixed integer/double data
------------------------------------------------------------------------- */

void CommLattice::perform_swap_general(Swap *swap)
{
  int i,j,m,n;
  int *index,*ivector;
  double *buf,*dvector;

  // post receives

  int ntotal = ninteger + ndouble;
  for (i = 0; i < swap->nrecv; i++)
    MPI_Irecv(swap->rdbuf[i],ntotal*swap->rcount[i],MPI_DOUBLE,
	      swap->rproc[i],0,world,&swap->request[i]);

  // pack data to send to each proc and send it

  for (i = 0; i < swap->nsend; i++) {
    index = swap->sindex[i];
    buf = swap->sdbuf;
    m = 0;
    for (n = 0; n < ninteger; n++) {
      ivector = iarray[n];
      for (j = 0; j < swap->scount[i]; j++) buf[m++] = ivector[index[j]];
    }
    for (n = 0; n < ndouble; n++) {
      dvector = darray[n];
      for (j = 0; j < swap->scount[i]; j++) buf[m++] = dvector[index[j]];
    }
    MPI_Send(buf,ntotal*swap->scount[i],MPI_DOUBLE,swap->sproc[i],0,world);
  }

  // wait on incoming messages

  if (swap->nrecv) MPI_Waitall(swap->nrecv,swap->request,swap->status);

  // unpack received buffers of data from each proc

  for (i = 0; i < swap->nrecv; i++) {
    index = swap->rindex[i];
    buf = swap->rdbuf[i];
    m = 0;
    for (n = 0; n < ninteger; n++) {
      ivector = iarray[n];
      for (j = 0; j < swap->rcount[i]; j++)
	ivector[index[j]] = static_cast<int> (buf[m++]);
    }
    for (n = 0; n < ndouble; n++) {
      dvector = darray[n];
      for (j = 0; j < swap->rcount[i]; j++) dvector[index[j]] = buf[m++];
    }
  }
}
