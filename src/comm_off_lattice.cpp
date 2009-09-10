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

#include "comm_off_lattice.h"
#include "app_off_lattice.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#include <map>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define DELTA 16384

enum{INTERIOR,EDGE,GHOST};        // same as in app_off_lattice.cpp

/* ---------------------------------------------------------------------- */

CommOffLattice::CommOffLattice(SPPARKS *spk) : Pointers(spk)
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

CommOffLattice::~CommOffLattice()
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
   if sectorflag = 0, just ghosts for entire proc domain 
   if sectorflag = 1, do all and sector ghosts, reverse if needed
   array = NULL = communicate iarray/darray from app
   array = non-NULL = communicate passed-in array (from diagnostic)
------------------------------------------------------------------------- */

void CommOffLattice::init(int nsector_request)
{
  appoff = (AppOffLattice *) app;

  ninteger = appoff->ninteger;
  ndouble = appoff->ndouble;

  xprd = appoff->xprd;
  yprd = appoff->yprd;
  zprd = appoff->zprd;

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
  //reverseswap = create_swap_all_reverse();

  /*
  nsector = nsector_request;
  if (nsector > 1) {
    sectorswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorswap[i] = create_swap_sector(appoff->set[i].nlocal,
					 appoff->set[i].site2i);
  }
  if (delreverse && nsector > 1) {
    sectorreverseswap = new Swap*[nsector];
    for (int i = 0; i < nsector; i++)
      sectorreverseswap[i] = 
	create_swap_sector_reverse(appoff->set[i].nlocal,
				   appoff->set[i].site2i);
  }
  */
}

/* ----------------------------------------------------------------------
   acquire ghost values for entire proc sub-domain
------------------------------------------------------------------------- */

void CommOffLattice::all()
{
  perform_swap(allswap);
}

/* ----------------------------------------------------------------------
   reverse communicate changed border values for entire proc sub-domain
------------------------------------------------------------------------- */

void CommOffLattice::all_reverse()
{
  //if (delreverse == 0.0) return;
  perform_swap(reverseswap);
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
  perform_swap(sectorreverseswap[isector]);
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost bins for my entire subdomain
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_all()
{
  int i,m;

  Swap *swap = new Swap;

  swap->nsend = 0;
  swap->sproc = NULL;
  swap->scount = NULL;
  swap->sindex = NULL;

  swap->nrecv = 0;
  swap->rproc = NULL;
  swap->rcount = NULL;
  swap->rindex = NULL;

  swap->request = NULL;
  swap->status = NULL;

  // self copy bins across PBC

  swap->ncopy = 0;
  swap->ccount = 0;
  swap->cbinsrc = NULL;
  swap->cbindest = NULL;

  int nbins = appoff->nbins;
  int *binflag = appoff->binflag;
  int *nbinimages = appoff->nbinimages;
  int **binimage = appoff->binimage;

  int ccount = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE) ccount += nbinimages[m];

  if (ccount) {
    swap->ncopy = 1;
    swap->ccount = ccount;
    swap->cbinsrc = new int[ccount];
    swap->cbindest = new int[ccount];
  }

  ccount = 0;
  for (m = 0; m < nbins; m++)
    if (binflag[m] == EDGE) {
      for (i = 0; i < nbinimages[m]; i++) {
	swap->cbinsrc[ccount] = m;
	swap->cbindest[ccount] = binimage[m][i];
	ccount++;
      }
    }

  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   reverse comm for my entire subdomain
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_all_reverse()
{
  Swap *swap = new Swap;
  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   acquire ghost bins for a single sector of my subdomain
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_sector(int nsites, 
							 int *site2i)
{
  Swap *swap = new Swap;
  return swap;
}

/* ----------------------------------------------------------------------
   create a Swap communication pattern
   reverse comm of a single sector
------------------------------------------------------------------------- */

CommOffLattice::Swap *CommOffLattice::create_swap_sector_reverse(int nsites,
								 int *site2i)
{
  Swap *swap = new Swap;
  return swap;
}

/* ----------------------------------------------------------------------
   free a Swap object
------------------------------------------------------------------------- */

void CommOffLattice::free_swap(Swap *swap)
{
  delete [] swap->sproc;
  delete [] swap->scount;
  for (int i = 0; i < swap->nsend; i++) memory->sfree(swap->sindex[i]);
  delete [] swap->sindex;

  delete [] swap->rproc;
  delete [] swap->rcount;
  for (int i = 0; i < swap->nrecv; i++) memory->sfree(swap->rindex[i]);
  delete [] swap->rindex;

  delete [] swap->cbinsrc;
  delete [] swap->cbindest;

  delete [] swap->request;
  delete [] swap->status;

  delete swap;
}

/* ----------------------------------------------------------------------
   communicate site values via Swap instructions
   setup per-site arrays for new ghost sites: bin,next,prev,nextimage
------------------------------------------------------------------------- */

void CommOffLattice::perform_swap(Swap *swap)
{
  int i,j,m,ibin,jbin,oldmax;

  // acquire bin info

  int *binhead = appoff->binhead;
  int **binoffset = appoff->binoffset;

  // unset nextimage() for all owned sites

  int *nextimage = appoff->nextimage;
  int nlocal = appoff->nlocal;
  for (int i = 0; i < nlocal; i++) nextimage[i] = -1;

  // delete all ghosts

  appoff->delete_all_ghosts();

  // nghost = # of ghosts I will receive
  // communicate # of sites in each bin I will send/receive
  // count # of sites in each on-processor copy bin

  int *next = appoff->next;
  int nghost = 0;

  if (swap->ncopy)
    for (m = 0; m < swap->ccount; m++) {
      i = binhead[swap->cbinsrc[m]];
      while (i >= 0) {
	nghost++;
	i = next[i];
      }
    }

  // grow per-site arrays as needed

  int nmax = appoff->nlocal + nghost;
  while (nmax > appoff->nmax) {
    oldmax = appoff->nmax;
    appoff->grow(0);
    appoff->add_free(oldmax);
  }

  // get current per-site pointers

  int *id = appoff->id;
  double **xyz = appoff->xyz;
  int *bin = appoff->bin;

  int *site = appoff->site;

  next = appoff->next;
  nextimage = appoff->nextimage;

  // post receives

  // perform on-processor copies across PBCs
  // add ghost sites one at a time
  // setup nextimage chain of ptrs for these sites

  if (swap->ncopy)
    for (m = 0; m < swap->ccount; m++) {
      ibin = swap->cbinsrc[m];
      jbin = swap->cbindest[m];
      i = binhead[ibin];
      while (i >= 0) {
	j = appoff->new_ghost();
	id[j] = id[i];
	xyz[j][0] = xyz[i][0] + binoffset[jbin][0]*xprd;
	xyz[j][1] = xyz[i][1] + binoffset[jbin][1]*yprd;
	xyz[j][2] = xyz[i][2] + binoffset[jbin][2]*zprd;
	site[j] = site[i];
	bin[j] = jbin;
	appoff->add_to_bin(j,jbin);
	nextimage[j] = nextimage[i];
	nextimage[i] = j;
	i = next[i];
      }
    }

  // wait on incoming messages

  // unpack received buffers of data from each proc
}
