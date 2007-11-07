/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice.h"
#include "app_lattice.h"
#include "sweep_lattice.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

#include <map>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 100

/* ---------------------------------------------------------------------- */

CommLattice::CommLattice(class SPK *spk) : SysPtr(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  allswap = NULL;
  sectorswap = NULL;
  nsector = 0;
}

/* ---------------------------------------------------------------------- */

CommLattice::~CommLattice()
{
  if (allswap) free_swap(allswap);
  for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
  delete [] sectorswap;
}

/* ---------------------------------------------------------------------- */

void CommLattice::init(SweepLattice *sweep,
		       const int delghost_in, const int dellocal_in)
{
  delghost = delghost_in;
  dellocal = dellocal_in;

  if (delghost > 1)
    error->all("More than 1st neighbor comm not yet supported by AppLattice");
  if (dellocal)
    error->all("Reverse comm not yet supported by AppLattice");

  AppLattice *applattice = (AppLattice *) app;
  dimension = applattice->dimension;
  int nlocal = applattice->nlocal;
  int *id = applattice->id;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  allswap = create_swap(nlocal,NULL,nlocal,id,numneigh,neighbor);

  if (sweep) {
    if (dimension == 2) nsector = 4;
    else nsector = 8;
    sectorswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorswap[i] = create_swap(nlocal,NULL,nlocal,id,numneigh,neighbor);
  }
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector
------------------------------------------------------------------------- */

void CommLattice::sector(int *lattice, const int isector)
{
  perform_swap(sectorswap[isector],lattice);
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
------------------------------------------------------------------------- */

void CommLattice::all(int *lattice)
{
  perform_swap(allswap,lattice);
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
------------------------------------------------------------------------- */

CommLattice::Swap *CommLattice::create_swap(int nsites, int *site2i,
					    int nlocal, int *id,
					    int *numneigh, int **neighbor)
{
  int i,j,k,m,size,original_proc;

  Swap *swap = new Swap;

  // loop over all neighbors of sites in site2i list
  // if site2i is NULL, then nsites = nlocal and loop over all owned sites
  // if I don't own neighbor and not already in ghost list,
  // add it to ghost list and to map

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;
  Ghost *buf = NULL;
  int maxghost = 0;
  int nghost = 0;

  for (m = 0; m < nsites; m++) {
    if (site2i) i = site2i[m];
    else i = m;
    for (j = 0; j < numneigh[i]; j++) {
      k = neighbor[i][j];
      if (k < nlocal) continue;
      if (hash.find(k) == hash.end()) {
	if (nghost == maxghost) {
	  maxghost += DELTA;
	  buf = (Ghost *) 
	    memory->srealloc(buf,maxghost*sizeof(Ghost),"comm:buf");
	}
	buf[nghost].id = id[k];
	buf[nghost].index = k;
	buf[nghost].proc = -1;
	hash.insert(std::pair<int,int> (k,0));
	nghost++;
      }
    }
  }

  // clear hash and put my entire list of owned lattice IDs in it

  hash.clear();
  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // maxghost = max ghosts on any proc

  int maxsize;
  MPI_Allreduce(&nghost,&maxsize,1,MPI_INT,MPI_MAX,world);

  buf = (Ghost *) memory->srealloc(buf,maxsize*sizeof(Ghost),"comm:buf");
  Ghost *bufcopy = (Ghost *) 
    memory->smalloc(maxsize*sizeof(Ghost),"comm:bufcopy");

  // initialize send lists

  int nsend = 0;
  int *sproc = NULL;
  int *scount = NULL;
  int **sindex = NULL;

  // cycle ghost list around ring of procs back to self
  // when receive it, fill in proc of sites I own
  // also add the ghost to my send lists

  MPI_Request request;
  MPI_Status status;

  size = nghost;
  original_proc = me;

  for (int loop = 0; loop < nprocs; loop++) {
    original_proc--;
    if (original_proc < 0) original_proc = nprocs-1;

    if (me != next) {
      MPI_Irecv(bufcopy,maxsize*sizeof(Ghost),MPI_CHAR,prev,0,world,&request);
      MPI_Send(buf,size*sizeof(Ghost),MPI_CHAR,next,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&size);
      size /= sizeof(Ghost);
      memcpy(buf,bufcopy,size*sizeof(Ghost));
    }
    for (i = 0; i < size; i++) {
      if (buf[i].proc >= 0) continue;
      loc = hash.find(buf[i].id);
      if (loc != hash.end()) {
	buf[i].proc = me;
	// add to send list based on j = loc->second;
	// if original_proc not on top of sproc, add it
	// add to sindex even if already there
	// need to avoid sending it multiple times to same proc
      }
    }
  }

  // original ghost list came back to me around ring
  // extract info for recv lists
  // error if any site is not filled in

  int nrecv = 0;
  int *rproc = NULL;
  int *rcount = NULL;
  int **rindex = NULL;

  for (i = 0; i < nghost; i++) {
    if (buf[i].proc == -1) error->one("Ghost site was not found");
  }

  // fill in swap data struct

  swap->nsend = nsend;
  swap->nrecv = nrecv;

  swap->sproc = sproc;
  swap->scount = scount;
  swap->sindex = sindex;
  // figure out sbuf size
  swap->sbuf = NULL;

  swap->rproc = rproc;
  swap->rcount = rcount;
  swap->rindex = rindex;
  // allocate rbufs
  swap->rbuf = NULL;

  if (swap->nrecv) {
    swap->request = new MPI_Request[swap->nrecv];
    swap->status = new MPI_Status[swap->nrecv];
  } else {
    swap->request = NULL;
    swap->status = NULL;
  }

  return swap;
}

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
------------------------------------------------------------------------- */

void CommLattice::free_swap(Swap *swap)
{
  delete [] swap->sproc;
  delete [] swap->scount;
  for (int i = 0; i < swap->nsend; i++) delete [] swap->sindex[i];
  delete [] swap->sindex;
  delete [] swap->sbuf;

  delete [] swap->rproc;
  delete [] swap->rcount;
  for (int i = 0; i < swap->nrecv; i++) delete [] swap->rindex[i];
  delete [] swap->rindex;
  for (int i = 0; i < swap->nrecv; i++) delete [] swap->rbuf[i];
  delete [] swap->rbuf;

  delete [] swap->request;
  delete [] swap->status;
}

/* ----------------------------------------------------------------------
   communicate ghost values via Swap instructions
------------------------------------------------------------------------- */

void CommLattice::perform_swap(Swap *swap, int *lattice)
{
  int i,j;
  int *index;
  int *buf;

  // post receives

  for (i = 0; i < swap->nrecv; i++)
    MPI_Irecv(swap->rbuf[i],swap->rcount[i],MPI_INT,swap->rproc[i],0,world,
	      &swap->request[i]);

  // pack data to send to each proc and send it

  for (i = 0; i < swap->nsend; i++) {
    index = swap->sindex[i];
    buf = swap->sbuf;
    for (j = 0; j < swap->scount[i]; j++) buf[j] = lattice[index[j]];
    MPI_Send(buf,swap->scount[i],MPI_INT,swap->sproc[i],0,world);
  }

  // wait on incoming messages

  if (swap->nrecv) MPI_Waitall(swap->nrecv,swap->request,swap->status);

  // unpack received buffers of data from each proc

  for (i = 0; i < swap->nrecv; i++) {
    index = swap->rindex[i];
    buf = swap->rbuf[i];
    for (j = 0; j < swap->rcount[i]; j++) lattice[index[j]] = buf[j];
  }
}
