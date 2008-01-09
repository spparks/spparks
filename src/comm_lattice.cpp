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

CommLattice::CommLattice(SPK *spk) : SysPtr(spk)
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
  int dimension = applattice->dimension;
  int nlocal = applattice->nlocal;
  int *id = applattice->id;
  int *numneigh = applattice->numneigh;
  int **neighbor = applattice->neighbor;

  if (allswap) free_swap(allswap);
  for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
  delete [] sectorswap;

  allswap = create_swap(nlocal,NULL,nlocal,id,numneigh,neighbor);

  if (sweep) {
    if (dimension == 2) nsector = 4;
    else nsector = 8;
    sectorswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorswap[i] = create_swap(sweep->sector[i].nlocal,
				  sweep->sector[i].site2i,
				  nlocal,id,numneigh,neighbor);
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
  int i,j,k,m,size,original_proc,isend,irecv;

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
	buf[nghost].id_global = id[k];
	buf[nghost].index_local = k;
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

  // cycle ghost list around ring of procs back to self
  // when receive it, fill in proc field for sites I own
  // also add this object to my send lists
  // original_proc = proc to send info to during swap

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
      loc = hash.find(buf[i].id_global);
      if (loc != hash.end()) {
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
	  sindex[isend] =
	    (int *) memory->srealloc(sindex[isend],smax[isend]*sizeof(int),
				     "comm:sindex");
	}
	sindex[isend][scount[isend]] = loc->second;
	scount[isend]++;
      }
    }
  }

  // clear hash

  hash.clear();

  // allocate recv lists

  int nrecv = 0;
  int *rproc = new int[nprocs];
  int *rcount = new int[nprocs];
  int *rmax = new int[nprocs];
  int **rindex = new int*[nprocs];
  int **rbuf = new int*[nprocs];
  for (int i = 0; i < nprocs; i++) {
    rmax[i] = 0;
    rindex[i] = NULL;
    rbuf[i] = NULL;
  }

  // original ghost list came back to me around ring
  // extract info for recv lists
  // error if any site is not filled in

  for (i = 0; i < nghost; i++) {
    if (buf[i].proc == -1) error->one("Ghost site was not found");

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
    
    // add index to send list going to a particular proc
    // grow rindex array if necessary

    if (rcount[irecv] == rmax[irecv]) {
      rmax[irecv] += DELTA;
      rindex[irecv] =
	(int *) memory->srealloc(rindex[irecv],rmax[irecv]*sizeof(int),
				 "comm:rindex");
    }
    rindex[irecv][rcount[irecv]] = buf[i].index_local;
    rcount[irecv]++;
  }

  // fill in swap data struct
  // allocate sbuf and rbuf[i] directly

  swap->nsend = nsend;
  swap->nrecv = nrecv;

  swap->sproc = sproc;
  swap->scount = scount;
  swap->smax = smax;
  swap->sindex = sindex;

  int max = 0;
  for (i = 0; i < nsend; i++) max = MAX(max,scount[i]);
  if (max) swap->sbuf = new int[max];
  else swap->sbuf = NULL;

  swap->rproc = rproc;
  swap->rcount = rcount;
  swap->rmax = rmax;
  swap->rindex = rindex;
  swap->rbuf = rbuf;

  for (i = 0; i < nrecv; i++)
    swap->rbuf[i] = new int[rcount[i]];

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
  delete [] swap->smax;
  for (int i = 0; i < swap->nsend; i++) memory->sfree(swap->sindex[i]);
  delete [] swap->sindex;
  delete [] swap->sbuf;

  delete [] swap->rproc;
  delete [] swap->rcount;
  delete [] swap->rmax;
  for (int i = 0; i < swap->nrecv; i++) memory->sfree(swap->rindex[i]);
  delete [] swap->rindex;
  for (int i = 0; i < swap->nrecv; i++) delete [] swap->rbuf[i];
  delete [] swap->rbuf;

  delete [] swap->request;
  delete [] swap->status;

  delete swap;
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
