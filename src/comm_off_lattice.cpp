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
#include "comm_off_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define BUFFACTOR 1.5

enum{INTERIOR,EDGE,GHOST};        // same as in app_off_lattice.cpp

/* ---------------------------------------------------------------------- */

CommOffLattice::CommOffLattice(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  allswap = NULL;
  sectorswap = NULL;
  sectorreverseswap = NULL;
  nsector = 0;

  smax = rmax = 0;
  sbuf = rbuf = NULL;
}

/* ---------------------------------------------------------------------- */

CommOffLattice::~CommOffLattice()
{
  if (allswap) free_swap(allswap);
  if (sectorswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
    delete [] sectorswap;
  }
  if (sectorreverseswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorreverseswap[i]);
    delete [] sectorreverseswap;
  }

  memory->destroy(sbuf);
  memory->destroy(rbuf);
}

/* ----------------------------------------------------------------------
   setup comm pattern
   3 kinds of patterns: all ghosts, sector ghosts, sector reverse
------------------------------------------------------------------------- */

void CommOffLattice::init(int nsector_request)
{
  appoff = (AppOffLattice *) app;

  // size_one = # of quantities per site = id + xyz + ninteger + ndouble

  ninteger = app->ninteger;
  ndouble = app->ndouble;
  size_one = 4 + ninteger + ndouble;
  if (ninteger == 1 && ndouble == 0) site_only = 1;
  else site_only = 0;

  xprd = domain->xprd;
  yprd = domain->yprd;
  zprd = domain->zprd;

  // clear out old swaps

  if (allswap) free_swap(allswap);
  if (sectorswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorswap[i]);
    delete [] sectorswap;
  }
  if (sectorreverseswap) {
    for (int i = 0; i < nsector; i++) free_swap(sectorreverseswap[i]);
    delete [] sectorreverseswap;
  }
  allswap = NULL;
  sectorswap = NULL;
  sectorreverseswap = NULL;

  // create new swaps as requested

  allswap = create_swap_all();

  nsector = nsector_request;
  if (nsector > 1) {
    setup_sector_chunks(appoff->nbinx,appoff->nbiny,appoff->nbinz);
    sectorswap = new Swap*[nsector];
    sectorreverseswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++) {
      sectorswap[i] = create_swap_sector(i);
      sectorreverseswap[i] = create_swap_sector_reverse(i);
    }
  }
}

/* ----------------------------------------------------------------------
   acquire ghost values for entire proc sub-domain
------------------------------------------------------------------------- */

void CommOffLattice::all()
{
  perform_swap(allswap);
}

/* ----------------------------------------------------------------------
   acquire ghost values for one sector
------------------------------------------------------------------------- */

void CommOffLattice::sector(int isector)
{
  perform_swap(sectorswap[isector]);
}

/* ----------------------------------------------------------------------
   reverse communicate changed border values for one sector
------------------------------------------------------------------------- */

void CommOffLattice::reverse_sector(int isector)
{
  perform_swap_reverse(sectorreverseswap[isector]);
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost bins for my entire subdomain
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_all()
{
  int i,m,n;

  int nbins = appoff->nbins;
  int *binflag = appoff->binflag;
  int *nimages = appoff->nimages;
  int **imageindex = appoff->imageindex;
  int **imageproc = appoff->imageproc;

  // add all image bins to list

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE) n += nimages[m];

  // allocate list

  int **list;
  memory->create(list,n,3,"commoff:list");

  // fill list with info on all image bins

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE)
      for (i = 0; i < nimages[m]; i++) {
	list[n][0] = m;
	list[n][1] = imageindex[m][i];
	list[n][2] = imageproc[m][i];
	n++;
      } 

  // create swap based on list of bins to send

  Swap *swap = new Swap;

  create_send_from_list(n,list,swap);
  create_recv_from_send(swap);

  memory->destroy(list);

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost bins for a single sector of my subdomain
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_sector(int isector)
{
  int i,m,n;

  int nbins = appoff->nbins;
  int nbinx = appoff->nbinx;
  int nbiny = appoff->nbiny;
  int nbinz = appoff->nbinz;
  int *nimages = appoff->nimages;
  int *binflag = appoff->binflag;
  int **imageindex = appoff->imageindex;
  int **imageproc = appoff->imageproc;

  // add only image bins to list needed for isector

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE)
      for (i = 0; i < nimages[m]; i++)
	if (bin_sector_ghost(isector,imageindex[m][i],nbinx,nbiny,nbinz)) n++;

  // allocate list

  int **list;
  memory->create(list,n,3,"commoff:list");

  // fill list with info on all image bins for isector

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE)
      for (i = 0; i < nimages[m]; i++)
	if (bin_sector_ghost(isector,imageindex[m][i],nbinx,nbiny,nbinz)) {
	  list[n][0] = m;
	  list[n][1] = imageindex[m][i];
	  list[n][2] = imageproc[m][i];
	  n++;
	}

  // create swap based on list of bins to send

  Swap *swap = new Swap;

  create_send_from_list(n,list,swap);
  create_recv_from_send(swap);

  memory->destroy(list);

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   reverse comm of a single sector
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_sector_reverse(int isector)
{
  int m,n;

  int nbins = appoff->nbins;
  int nbinx = appoff->nbinx;
  int nbiny = appoff->nbiny;
  int nbinz = appoff->nbinz;
  int *binflag = appoff->binflag;
  int *ghostindex = appoff->ghostindex;
  int *ghostproc = appoff->ghostproc;

  // add only ghost bins to list needed for isector

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == GHOST)
      if (bin_sector_ghost(isector,m,nbinx,nbiny,nbinz)) n++;

  // allocate list

  int **list;
  memory->create(list,n,3,"commoff:list");

  // fill list with info on all ghost bins for isector

  n = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == GHOST)
      if (bin_sector_ghost(isector,m,nbinx,nbiny,nbinz)) {
	list[n][0] = m;
	list[n][1] = ghostindex[m];
	list[n][2] = ghostproc[m];
	n++;
      }

  // create swap based on list of bins to send

  Swap *swap = new Swap;

  create_send_from_list(n,list,swap);
  create_recv_from_send(swap);

  memory->destroy(list);

  return swap;
}

/* ----------------------------------------------------------------------
   create send and copy portions of a Swap communication pattern
   generated from list of bins I send to others or self
------------------------------------------------------------------------- */

void CommOffLattice::create_send_from_list(int nlist, int **list, Swap *swap)
{
  int i,j,isend;

  // proc[i] = count of bins I send to proc I, including self

  int *proc = new int[nprocs];
  for (i = 0; i < nprocs; i++) proc[i] = 0;
  for (i = 0; i < nlist; i++) proc[list[i][2]]++;

  // nsend = # of procs I send to, excluding self

  int nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (proc[i]) nsend++;
  if (proc[me]) nsend--;

  // ccount = # of bins I copy to myself
  // cbinsrc,cbindest = list of bins copied from src to dest

  int ccount = proc[me];
  int *cbinsrc = NULL;
  int *cbindest = NULL;
  if (ccount) {
    cbinsrc = new int[ccount];
    cbindest = new int[ccount];
  }

  ccount = 0;
  for (i = 0; i < nlist; i++) {
    if (list[i][2] == me) {
      cbinsrc[ccount] = list[i][0];
      cbindest[ccount] = list[i][1];
      ccount++;
    }
  }

  // setup send lists
  // sproc = list of procs I send to, excluding self
  // scount = # of bins I send to each proc
  // start sproc with procs beyond me, so comm is more staggered

  int *sproc = new int[nsend];
  int *scount = new int[nsend];
  int *stotal = new int[nsend];

  i = me;
  nsend = 0;
  for (j = 0; j < nprocs-1; j++) {
    i++;
    if (i == nprocs) i = 0;
    if (proc[i]) {
      sproc[nsend] = i;
      scount[nsend] = proc[i];
      proc[i] = nsend;
      nsend++;
    }
  }

  // allocate and setup ragged arrays: sindex and ssite
  // sindex = list of bins to send
  // ssite = messages to send with IDs of bins for receiving procs
  // ssite will actually be sent in create_recv()
  // use proc[] to convert proc ID to isend

  int **sindex;
  memory->create_ragged(sindex,nsend,scount,"commoff:sindex");
  int **ssite;
  memory->create_ragged(ssite,nsend,scount,"commoff:ssite");

  for (isend = 0; isend < nsend; isend++) scount[isend] = 0;

  for (i = 0; i < nlist; i++) {
    if (list[i][2] != me) {
      isend = proc[list[i][2]];
      sindex[isend][scount[isend]] = list[i][0];
      ssite[isend][scount[isend]] = list[i][1];
      scount[isend]++;
    }
  }

  // clean up

  delete [] proc;

  // fill in swap data structure with send and copy info

  swap->nsend = nsend;
  swap->sproc = sproc;
  swap->scount = scount;
  swap->sindex = sindex;
  swap->ssite = ssite;
  swap->stotal = stotal;

  if (ccount) swap->ncopy = 1;
  else swap->ncopy = 0;
  swap->ccount = ccount;
  swap->cbinsrc = cbinsrc;
  swap->cbindest = cbindest;
}

/* ----------------------------------------------------------------------
   create recv portion of a Swap communication pattern
   communicate site list from sending procs to populate recv lists
------------------------------------------------------------------------- */

void CommOffLattice::create_recv_from_send(Swap *swap)
{
  int i,isend;

  // swap params

  int nsend = swap->nsend;
  int *sproc = swap->sproc;
  int *scount = swap->scount;
  int **ssite = swap->ssite;

  // nrecv = # of procs I receive from

  int *proc = new int[nprocs];
  int *count = new int[nprocs];

  for (i = 0; i < nprocs; i++) {
    proc[i] = 0;
    count[i] = 1;
  }
  for (isend = 0; isend < nsend; isend++) proc[sproc[isend]] = 1;

  int nrecv;
  MPI_Reduce_scatter(proc,&nrecv,count,MPI_INT,MPI_SUM,world);

  // setup recv lists

  int *rproc = new int[nrecv];
  int *rcount = new int[nrecv];
  int *rtotal = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // rproc = list of procs I receive from
  // rcount = list of procs I receive from
  // rcount = # of bins I recv from each proc
  // fill these by telling receivers how many bins I send

  for (i = 0; i < nsend; i++)
    MPI_Send(&scount[i],1,MPI_INT,sproc[i],0,world);

  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&rcount[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    rproc[i] = status->MPI_SOURCE;
  }

  // sort rproc and rcount in ascending rproc order
  // this insures same answer from run to run (on same # of procs)
  // else random ordering of MPI_ANY_SOURCE could induce different
  //   ordering of my owned sites as reverse comm occurs

  int j,tmp;
  for (i = 0; i < nrecv; i++)
    for (j = i+1; j < nrecv; j++) {
      if (rproc[j] < rproc[i]) {
	tmp = rproc[j];
	rproc[j] = rproc[i];
	rproc[i] = tmp;
	tmp = rcount[j];
	rcount[j] = rcount[i];
	rcount[i] = tmp;
      }
    }

  MPI_Barrier(world);

  // allocate and setup ragged arrays: rindex and rsite
  // rindex = list of bins to receive
  // fill rindex by sending ssite to receivers

  int **rindex;
  memory->create_ragged(rindex,nrecv,rcount,"commoff:rindex");
  int **rsite;
  memory->create_ragged(rsite,nrecv,rcount,"commoff:rsite");

  for (int irecv = 0; irecv < nrecv; irecv++)
    MPI_Irecv(rindex[irecv],rcount[irecv],MPI_INT,
	      rproc[irecv],0,world,&request[irecv]);

  MPI_Barrier(world);

  for (int isend = 0; isend < nsend; isend++)
    MPI_Send(ssite[isend],scount[isend],MPI_INT,
	     sproc[isend],0,world);

  if (nrecv) MPI_Waitall(nrecv,request,status);

  // clean up

  delete [] proc;
  delete [] count;

  // fill in swap data structure with recv info

  swap->nrecv = nrecv;
  swap->rproc = rproc;
  swap->rcount = rcount;
  swap->rindex = rindex;
  swap->rsite = rsite;
  swap->rtotal = rtotal;
  swap->request = request;
  swap->status = status;
}

/* ----------------------------------------------------------------------
   free a Swap object
------------------------------------------------------------------------- */

void CommOffLattice::free_swap(Swap *swap)
{
  delete [] swap->sproc;
  delete [] swap->scount;
  memory->destroy(swap->sindex);
  memory->destroy(swap->ssite);
  delete [] swap->stotal;

  delete [] swap->rproc;
  delete [] swap->rcount;
  memory->destroy(swap->rindex);
  memory->destroy(swap->rsite);
  delete [] swap->rtotal;

  delete [] swap->cbinsrc;
  delete [] swap->cbindest;

  delete [] swap->request;
  delete [] swap->status;

  delete swap;
}

/* ----------------------------------------------------------------------
   communicate site values via Swap instructions
   sites from owned bins will be communicated or copied
   sites in ghost bins will be created
------------------------------------------------------------------------- */

void CommOffLattice::perform_swap(Swap *swap)
{
  int i,j,k,m,n,ibin,jbin,oldmax,count,total;
  int *site;

  // swap info

  int nsend = swap->nsend;
  int nrecv = swap->nrecv;
  int *sproc = swap->sproc;
  int *scount = swap->scount;
  int **sindex = swap->sindex;
  int **ssite = swap->ssite;
  int *stotal = swap->stotal;
  int *rproc = swap->rproc;
  int *rcount = swap->rcount;
  int **rindex = swap->rindex;
  int **rsite = swap->rsite;
  int *rtotal = swap->rtotal;
  int ncopy = swap->ncopy;
  int ccount = swap->ccount;
  int *cbinsrc = swap->cbinsrc;
  int *cbindest = swap->cbindest;
  MPI_Request *request = swap->request;
  MPI_Status *status = swap->status;

  // app info

  int *binhead = appoff->binhead;
  int **pbcoffset = appoff->pbcoffset;
  int *next = appoff->next;
  int *nextimage = appoff->nextimage;
  int nlocal = appoff->nlocal;

  // if no sectors on single proc, unset nextimage() for all owned sites

  if (nsector == 1)
    for (int i = 0; i < nlocal; i++) nextimage[i] = -1;

  // delete all ghosts, since are replacing them

  appoff->delete_all_ghosts();

  // post receives for bin counts

  for (int irecv = 0; irecv < nrecv; irecv++)
    MPI_Irecv(rsite[irecv],rcount[irecv],MPI_INT,
	      rproc[irecv],0,world,&request[irecv]);

  // barrier to insure receives are posted

  MPI_Barrier(world);

  // count sites in each bin I send
  // sendmax = max # of sites to send to one proc

  int sendmax = 0;

  for (int isend = 0; isend < nsend; isend++) {
    total = 0;
    for (n = 0; n < scount[isend]; n++) {
      count = 0;
      m = binhead[sindex[isend][n]];
      while (m >= 0) {
	count++;
	m = next[m];
      }
      ssite[isend][n] = count;
      total += count;
    }
    stotal[isend] = total;
    sendmax = MAX(sendmax,total);
  }

  // send bin counts

  for (int isend = 0; isend < nsend; isend++)
    MPI_Send(ssite[isend],scount[isend],MPI_INT,
	     sproc[isend],0,world);

  // copycount = # of ghosts generated via copying from my own bins

  int copycount = 0;

  if (ncopy)
    for (n = 0; n < ccount; n++) {
      m = binhead[cbinsrc[n]];
      while (m >= 0) {
	copycount++;
	m = next[m];
      }
    }

  // wait on all incoming bin counts

  if (nrecv) MPI_Waitall(nrecv,request,status);

  // recvcount = # of ghosts received from other procs

  int recvcount = 0;

  for (int irecv = 0; irecv < nrecv; irecv++) {
    total = 0;
    for (n = 0; n < rcount[irecv]; n++)
      total += rsite[irecv][n];
    rtotal[irecv] = total;
    recvcount += total;
  }

  // trigger app to grow per-site arrays if new ghosts exceed limit

  int nmax = nlocal + recvcount + copycount;
  while (nmax > appoff->nmax) {
    oldmax = appoff->nmax;
    appoff->grow(0);
    appoff->add_free(oldmax);
  }

  // reset site ptrs from app

  tagint *id = app->id;
  double **xyz = app->xyz;
  int **iarray = app->iarray;
  double **darray = app->darray;
  if (site_only) site = iarray[0];
  int *bin = appoff->bin;
  next = appoff->next;
  nextimage = appoff->nextimage;

  // grow send/receive buffers as needed

  if (sendmax > smax) {
    smax = static_cast<int> (BUFFACTOR*sendmax);
    memory->destroy(sbuf);
    memory->create(sbuf,smax*size_one,"app:sbuf");
  }

  if (recvcount > rmax) {
    rmax = static_cast<int> (BUFFACTOR*recvcount);
    memory->destroy(rbuf);
    memory->create(rbuf,rmax*size_one,"app:rbuf");
  }

  // post receives for sites

  int offset = 0;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    MPI_Irecv(&rbuf[offset],size_one*rtotal[irecv],MPI_DOUBLE,
	      rproc[irecv],0,world,&request[irecv]);
    offset += size_one*rtotal[irecv];
  }

  // barrier to insure receives are posted

  MPI_Barrier(world);

  // send sites

  for (int isend = 0; isend < nsend; isend++) {
    i = 0;
    for (n = 0; n < scount[isend]; n++) {
      m = binhead[sindex[isend][n]];
      while (m >= 0) {
	sbuf[i++] = id[m];
	sbuf[i++] = xyz[m][0];
	sbuf[i++] = xyz[m][1];
	sbuf[i++] = xyz[m][2];
	if (site_only) sbuf[i++] = site[m];
	else {
	  for (k = 0; k < ninteger; k++) sbuf[i++] = iarray[k][i];
	  for (k = 0; k < ndouble; k++) sbuf[i++] = darray[k][i];
	}
	m = next[m];
      }
    }

    MPI_Send(sbuf,size_one*stotal[isend],MPI_DOUBLE,sproc[isend],0,world);
  }

  // perform on-processor copies across PBCs
  // add ghost sites one at a time
  // if no sectors on single proc, setup nextimage chain of ptrs

  if (ncopy)
    for (m = 0; m < ccount; m++) {
      ibin = cbinsrc[m];
      jbin = cbindest[m];
      i = binhead[ibin];
      while (i >= 0) {
	j = appoff->new_ghost_site();
	id[j] = id[i];
	xyz[j][0] = xyz[i][0] + pbcoffset[jbin][0]*xprd;
	xyz[j][1] = xyz[i][1] + pbcoffset[jbin][1]*yprd;
	xyz[j][2] = xyz[i][2] + pbcoffset[jbin][2]*zprd;
	if (site_only) site[j] = site[i];
	else {
	  for (k = 0; k < ninteger; k++) iarray[k][j] = iarray[k][i];
	  for (k = 0; k < ndouble; k++) darray[k][j] = darray[k][i];
	}
	bin[j] = jbin;
	appoff->add_to_bin(j,jbin);
	if (nsector == 1) {
	  nextimage[j] = nextimage[i];
	  nextimage[i] = j;
	}
	i = next[i];
      }
    }

  // wait on all incoming sites

  if (nrecv) MPI_Waitall(nrecv,request,status);

  // unpack received sites from each proc
  // add sites to my bins, adding in PBC if necessary

  m = 0;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    for (n = 0; n < rcount[irecv]; n++) {
      for (i = 0; i < rsite[irecv][n]; i++) {
	j = appoff->new_ghost_site();
	jbin = rindex[irecv][n];
	id[j] = static_cast<tagint> (rbuf[m++]);
	xyz[j][0] = rbuf[m++] + pbcoffset[jbin][0]*xprd;
	xyz[j][1] = rbuf[m++] + pbcoffset[jbin][1]*yprd;
	xyz[j][2] = rbuf[m++] + pbcoffset[jbin][2]*zprd;
	if (site_only) site[j] = static_cast<int> (rbuf[m++]);
	else {
	  for (k = 0; k < ninteger; k++) 
	    iarray[k][j] = static_cast<int> (rbuf[m++]);
	  for (k = 0; k < ndouble; k++) darray[k][j] = rbuf[m++];
	}
	bin[j] = jbin;
	appoff->add_to_bin(j,jbin);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   communicate site values via Swap instructions
   owned sites int ghost bins will be communicated or copied
   sites in owned bins will be created
------------------------------------------------------------------------- */

void CommOffLattice::perform_swap_reverse(Swap *swap)
{
  int i,j,k,m,n,ibin,jbin,oldmax,count,total,nextptr;
  int *site;

  // swap info

  int nsend = swap->nsend;
  int nrecv = swap->nrecv;
  int *sproc = swap->sproc;
  int *scount = swap->scount;
  int **sindex = swap->sindex;
  int **ssite = swap->ssite;
  int *stotal = swap->stotal;
  int *rproc = swap->rproc;
  int *rcount = swap->rcount;
  int **rindex = swap->rindex;
  int **rsite = swap->rsite;
  int *rtotal = swap->rtotal;
  int ncopy = swap->ncopy;
  int ccount = swap->ccount;
  int *cbinsrc = swap->cbinsrc;
  int *cbindest = swap->cbindest;
  MPI_Request *request = swap->request;
  MPI_Status *status = swap->status;

  // app info

  int *binhead = appoff->binhead;
  int **pbcoffset = appoff->pbcoffset;
  int *next = appoff->next;
  int nlocal = appoff->nlocal;

  // post receives for bin counts

  for (int irecv = 0; irecv < nrecv; irecv++)
    MPI_Irecv(rsite[irecv],rcount[irecv],MPI_INT,
	      rproc[irecv],0,world,&request[irecv]);

  // barrier to insure receives are posted

  MPI_Barrier(world);

  // count owned sites in each ghost bin I send
  // sendmax = max # of sites to send to one proc
  // sendcount = total # of sites I send

  int sendmax = 0;
  int sendcount = 0;

  for (int isend = 0; isend < nsend; isend++) {
    total = 0;
    for (n = 0; n < scount[isend]; n++) {
      count = 0;
      m = binhead[sindex[isend][n]];
      while (m >= 0) {
	if (m < nlocal) count++;
	m = next[m];
      }
      ssite[isend][n] = count;
      total += count;
    }
    stotal[isend] = total;
    sendmax = MAX(sendmax,total);
    sendcount += total;
  }

  // send bin counts

  for (int isend = 0; isend < nsend; isend++)
    MPI_Send(ssite[isend],scount[isend],MPI_INT,
	     sproc[isend],0,world);

  // copycount = # of sites copied from my own bins

  int copycount = 0;

  if (ncopy)
    for (n = 0; n < ccount; n++) {
      m = binhead[cbinsrc[n]];
      while (m >= 0) {
	if (m < nlocal) copycount++;
	m = next[m];
      }
    }

  // wait on all incoming bin counts

  if (nrecv) MPI_Waitall(nrecv,request,status);

  // recvcount = # of owned sites received from other procs

  int recvcount = 0;

  for (int irecv = 0; irecv < nrecv; irecv++) {
    total = 0;
    for (n = 0; n < rcount[irecv]; n++)
      total += rsite[irecv][n];
    rtotal[irecv] = total;
    recvcount += total;
  }

  // trigger app to grow per-site arrays if new owned sites exceed limit

  int nmax = appoff->nlocal - sendcount + recvcount;
  while (nmax > appoff->nmax) {
    oldmax = appoff->nmax;
    appoff->grow(0);
    appoff->add_free(oldmax);
  }

  // reset site ptrs from app

  tagint *id = app->id;
  double **xyz = app->xyz;
  int **iarray = app->iarray;
  double **darray = app->darray;
  if (site_only) site = iarray[0];
  int *bin = appoff->bin;
  next = appoff->next;

  // grow send/receive buffers as needed

  if (sendmax > smax) {
    smax = static_cast<int> (BUFFACTOR*sendmax);
    memory->destroy(sbuf);
    memory->create(sbuf,smax*size_one,"app:sbuf");
  }

  if (recvcount > rmax) {
    rmax = static_cast<int> (BUFFACTOR*recvcount);
    memory->destroy(rbuf);
    memory->create(rbuf,rmax*size_one,"app:rbuf");
  }

  // post receives for sites

  int offset = 0;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    MPI_Irecv(&rbuf[offset],size_one*rtotal[irecv],MPI_DOUBLE,
	      rproc[irecv],0,world,&request[irecv]);
    offset += size_one*rtotal[irecv];
  }

  // barrier to insure receives are posted

  MPI_Barrier(world);

  // send sites, subtracting PBC if necessary
  // delete each owned site after packing it
  // update nlocal since delete_owned_site() will decrement it

  for (int isend = 0; isend < nsend; isend++) {
    i = 0;
    for (n = 0; n < scount[isend]; n++) {
      ibin = sindex[isend][n];
      m = binhead[ibin];
      while (m >= 0) {
	if (m < nlocal) {
	  sbuf[i++] = id[m];
	  sbuf[i++] = xyz[m][0] - pbcoffset[ibin][0]*xprd;
	  sbuf[i++] = xyz[m][1] - pbcoffset[ibin][1]*yprd;
	  sbuf[i++] = xyz[m][2] - pbcoffset[ibin][2]*zprd;
	  if (site_only) sbuf[i++] = site[m];
	  else {
	    for (k = 0; k < ninteger; k++) sbuf[i++] = iarray[k][m];
	    for (k = 0; k < ndouble; k++) sbuf[i++] = darray[k][m];
	  }
	  m = appoff->delete_owned_site(m);
	  nlocal = appoff->nlocal;
	} else m = next[m];
      }
    }

    MPI_Send(sbuf,size_one*stotal[isend],MPI_DOUBLE,sproc[isend],0,world);
  }

  // perform on-processor copies across PBC from ghost bin to edge bin

  if (ncopy)
    for (m = 0; m < ccount; m++) {
      ibin = cbinsrc[m];
      jbin = cbindest[m];
      i = binhead[ibin];
      while (i >= 0) {
	if (i < nlocal) {
	  nextptr = next[i];
	  xyz[i][0] -= pbcoffset[ibin][0]*xprd;
	  xyz[i][1] -= pbcoffset[ibin][1]*yprd;
	  xyz[i][2] -= pbcoffset[ibin][2]*zprd;
	  appoff->delete_from_bin(i,ibin);
	  appoff->add_to_bin(i,jbin);
	  bin[i] = jbin;
	  i = nextptr;
	} else i = next[i];
      }
    }

  // delete all ghosts, since have finished looping over ghost bin sites

  appoff->delete_all_ghosts();

  // wait on all incoming sites

  if (nrecv) MPI_Waitall(nrecv,request,status);

  // unpack received sites from each proc
  // add sites to my edge bins

  m = 0;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    for (n = 0; n < rcount[irecv]; n++) {
      for (i = 0; i < rsite[irecv][n]; i++) {
	j = appoff->new_owned_site();
	jbin = rindex[irecv][n];
	id[j] = static_cast<tagint> (rbuf[m++]);
	xyz[j][0] = rbuf[m++];
	xyz[j][1] = rbuf[m++];
	xyz[j][2] = rbuf[m++];
	if (site_only) site[j] = static_cast<int> (rbuf[m++]);
	else {
	  for (k = 0; k < ninteger; k++) 
	    iarray[k][j] = static_cast<int> (rbuf[m++]);
	  for (k = 0; k < ndouble; k++) darray[k][m] = rbuf[m++];
	}
	bin[j] = jbin;
	appoff->add_to_bin(j,jbin);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return 1 if ibin is a ghost of sector isector
   ibin is a local ghost bin ID
   nbinx,nbiny,nbinz = extent of local bins
------------------------------------------------------------------------- */

int CommOffLattice::bin_sector_ghost(int isector, int ibin, 
				     int nbinx, int nbiny, int nbinz)
{
  int i = ibin % nbinx;
  int j = (ibin/nbinx) % nbiny;
  int k = ibin / (nbinx*nbiny);

  for (int m = 0; m < nchunk; m++)
    if (i >= chunklo[isector][m][0] && i <= chunkhi[isector][m][0] &&
	j >= chunklo[isector][m][1] && j <= chunkhi[isector][m][1] &&
	k >= chunklo[isector][m][2] && k <= chunkhi[isector][m][2]) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   setup sector chunks
   nchunk = 3 planes for dim = 3, nsector = 8
   nchunk = 4 planes for dim = 3, nsector = 4
   nchunk = 5 planes for dim = 3, nsector = 2
   nchunk = 2 lines for dim = 2, nsector = 4
   nchunk = 3 lines for dim = 2, nsector = 2
------------------------------------------------------------------------- */

void CommOffLattice::setup_sector_chunks(int nbinx, int nbiny, int nbinz)
{
  // nbinx,nbiny,nbinz are all even

  int ihalf = nbinx/2;
  int jhalf = nbiny/2;
  int khalf = nbinz/2;

  // fill in chunk boundaries for dimension = 2/3, nsector = 2/4/8

  int dimension = appoff->dimension;

  if (dimension == 2 && nsector == 4) {
    nchunk = 2;
    fill_chunk(0,0,0,ihalf,0,0,0,0);
    fill_chunk(0,1,0,0,0,jhalf,0,0);
    fill_chunk(1,0,ihalf-1,nbinx-1,0,0,0,0);
    fill_chunk(1,1,nbinx-1,nbinx-1,0,jhalf,0,0);
    fill_chunk(2,0,0,ihalf,nbiny-1,nbiny-1,0,0);
    fill_chunk(2,1,0,0,jhalf-1,nbiny-1,0,0);
    fill_chunk(3,0,ihalf-1,nbinx-1,nbiny-1,nbiny-1,0,0);
    fill_chunk(3,1,nbinx-1,nbinx-1,jhalf-1,nbiny-1,0,0);
  }

  if (dimension == 2 && nsector == 2) {
    nchunk = 3;
    fill_chunk(0,0,0,0,0,nbiny-1,0,0);
    fill_chunk(0,1,0,ihalf,0,0,0,0);
    fill_chunk(0,2,0,ihalf,nbiny-1,nbiny-1,0,0);
    fill_chunk(1,0,nbinx-1,nbinx-1,0,nbiny-1,0,0);
    fill_chunk(1,1,ihalf-1,nbinx-1,0,0,0,0);
    fill_chunk(1,2,ihalf-1,nbinx-1,nbiny-1,nbiny-1,0,0);
  }

  if (dimension == 3 && nsector == 8) {
    nchunk = 3;
    fill_chunk(0,0,0,ihalf,0,jhalf,0,0);
    fill_chunk(0,1,0,0,0,jhalf,0,khalf);
    fill_chunk(0,2,0,ihalf,0,0,0,khalf);
    fill_chunk(1,0,ihalf-1,nbinx-1,0,jhalf,0,0);
    fill_chunk(1,1,nbinx-1,nbinx-1,0,jhalf,0,khalf);
    fill_chunk(1,2,ihalf-1,nbinx-1,0,0,0,khalf);
    fill_chunk(2,0,0,ihalf,jhalf-1,nbiny-1,0,0);
    fill_chunk(2,1,0,0,jhalf-1,nbiny-1,0,khalf);
    fill_chunk(2,2,0,ihalf,nbiny-1,nbiny-1,0,khalf);
    fill_chunk(3,0,ihalf-1,nbinx-1,jhalf-1,nbiny-1,0,0);
    fill_chunk(3,1,nbinx-1,nbinx-1,jhalf-1,nbiny-1,0,khalf);
    fill_chunk(3,2,ihalf-1,nbinx-1,nbiny-1,nbiny-1,0,khalf);
    fill_chunk(4,0,0,ihalf,0,jhalf,nbinx-1,nbinz-1);
    fill_chunk(4,1,0,0,0,jhalf,khalf-1,nbinz-1);
    fill_chunk(4,2,0,ihalf,0,0,khalf-1,nbinz-1);
    fill_chunk(5,0,ihalf-1,nbinx-1,0,jhalf,nbinz-1,nbinz-1);
    fill_chunk(5,1,nbinx-1,nbinx-1,0,jhalf,khalf-1,nbinz-1);
    fill_chunk(5,2,ihalf-1,nbinx-1,0,0,khalf-1,nbinz-1);
    fill_chunk(6,0,0,ihalf,jhalf-1,nbiny-1,nbinz-1,nbinz-1);
    fill_chunk(6,1,0,0,jhalf-1,nbiny-1,khalf-1,nbinz-1);
    fill_chunk(6,2,0,ihalf,nbiny-1,nbiny-1,khalf-1,nbinz-1);
    fill_chunk(7,0,ihalf-1,nbinx-1,jhalf-1,nbiny-1,nbinz-1,nbinz-1);
    fill_chunk(7,1,nbinx-1,nbinx-1,jhalf-1,nbiny-1,khalf-1,nbinz-1);
    fill_chunk(7,2,ihalf-1,nbinx-1,nbiny-1,nbiny-1,khalf-1,nbinz-1);
  }

  if (dimension == 3 && nsector == 4) {
    nchunk = 4;
    fill_chunk(0,0,0,ihalf,0,0,0,nbinz-1);
    fill_chunk(0,1,0,0,0,jhalf,0,nbinz-1);
    fill_chunk(0,2,0,ihalf,0,jhalf,0,0);
    fill_chunk(0,3,0,ihalf,0,jhalf,nbinz-1,nbinz-1);
    fill_chunk(1,0,ihalf-1,nbinx-1,0,0,0,nbinz-1);
    fill_chunk(1,1,nbinx-1,nbinx-1,0,jhalf,0,nbinz-1);
    fill_chunk(1,2,ihalf-1,nbinx-1,0,jhalf,0,0);
    fill_chunk(1,3,ihalf-1,nbinx-1,0,jhalf,nbinz-1,nbinz-1);
    fill_chunk(2,0,0,ihalf,nbiny-1,nbiny-1,0,nbinz-1);
    fill_chunk(2,1,0,0,jhalf-1,nbiny-1,0,nbinz-1);
    fill_chunk(2,2,0,ihalf,jhalf-1,nbiny-1,0,0);
    fill_chunk(2,3,0,ihalf,jhalf-1,nbiny-1,nbinz-1,nbinz-1);
    fill_chunk(3,0,ihalf-1,nbinx-1,nbiny-1,nbiny-1,0,nbinz-1);
    fill_chunk(3,1,nbinx-1,nbinx-1,jhalf-1,nbiny-1,0,nbinz-1);
    fill_chunk(3,2,ihalf-1,nbinx-1,jhalf-1,nbiny-1,0,0);
    fill_chunk(3,3,ihalf-1,nbinx-1,jhalf-1,nbiny-1,nbinz-1,nbinz-1);
  }

  if (dimension == 3 && nsector == 2) {
    nchunk = 5;
    fill_chunk(0,0,0,0,0,nbiny-1,0,nbinz-1);
    fill_chunk(0,1,0,ihalf,0,0,0,nbinz-1);
    fill_chunk(0,2,0,ihalf,nbiny-1,nbiny-1,0,nbinz-1);
    fill_chunk(0,3,0,ihalf,0,nbiny-1,0,0);
    fill_chunk(0,4,0,ihalf,0,nbiny-1,nbinz-1,nbinz-1);
    fill_chunk(1,0,nbinx-1,nbinx-1,0,nbiny-1,0,nbinz-1);
    fill_chunk(1,1,ihalf-1,nbinx-1,0,0,0,nbinz-1);
    fill_chunk(1,3,ihalf-1,nbinx-1,nbiny-1,nbiny-1,0,nbinz-1);
    fill_chunk(1,4,ihalf-1,nbinx-1,0,nbiny-1,0,0);
    fill_chunk(1,5,ihalf-1,nbinx-1,0,nbiny-1,nbinz-1,nbinz-1);
  }
}

/* ----------------------------------------------------------------------
   set bounds of a single chunk
------------------------------------------------------------------------- */

void CommOffLattice::fill_chunk(int isector, int ichunk,
				int ilo, int ihi, int jlo, int jhi, 
				int klo, int khi)
{
  chunklo[isector][ichunk][0] = ilo;
  chunklo[isector][ichunk][1] = jlo;
  chunklo[isector][ichunk][2] = klo;
  chunkhi[isector][ichunk][0] = ihi;
  chunkhi[isector][ichunk][1] = jhi;
  chunkhi[isector][ichunk][2] = khi;
}
