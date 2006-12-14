/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_grain.h"
#include "comm_grain.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

CommGrain::CommGrain(class SPK *spk) : SysPtr(spk)
{
  swapinfo = NULL;
}

/* ---------------------------------------------------------------------- */

CommGrain::~CommGrain()
{
  free_swap();
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommGrain::setup(const int nx, const int ny, 
		       const int nxh, const int nyh, 
		       const int procwest, const int proceast, 
		       const int procsouth, const int procnorth)
{
  // initialize swap parameters and allocate memory

  SwapInfo* swap;
  int iswap,isector;

  // isector = 0,1,2,3 refers to SW,NW,SE,NE quadrants, respectively
  // iswap = 0 refers to the north/south swap
  // iswap = 1 refers to the east/west swap 

  nsector = 4;
  nswap = 2;
  allocate_swap(nsector,nswap);

  // southwest sector; send north

  isector = 0;
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = 0;
  swap->recvnum = nxh+1;
  swap->sendix = 0;
  swap->sendiy = ny;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procsouth;
  swap->sendproc = procnorth;
  swap->stride = ny+2;
  swap->copy0ixa = nx;
  swap->copy0ixb = 0;

  // southwest sector; send east

  isector = 0;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = 0;
  swap->recvnum = nyh+1;
  swap->sendix = nx;
  swap->sendiy = 0;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procwest;
  swap->sendproc = proceast;
  swap->stride = 1;
  swap->copy0ixa = -1;
  swap->copy0ixb = -1;

  // northwest sector; send south

  isector = 1;
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = ny+1;
  swap->recvnum = nxh+1;
  swap->sendix = 0;
  swap->sendiy = 1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procnorth;
  swap->sendproc = procsouth;
  swap->stride = ny+2;
  swap->copy0ixa = nx;
  swap->copy0ixb = 0;

  // northwest sector; send east

  isector = 1;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = nyh-1;
  swap->recvnum = ny-nyh+3;
  swap->sendix = nx;
  swap->sendiy = nyh-1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procwest;
  swap->sendproc = proceast;
  swap->stride = 1;
  swap->copy0ixa = -1;
  swap->copy0ixb = -1;

  // southeast sector; send north

  isector = 2;
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nxh-1;
  swap->recviy = 0;
  swap->recvnum = nx-nxh+3;
  swap->sendix = nxh-1;
  swap->sendiy = ny;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procsouth;
  swap->sendproc = procnorth;
  swap->stride = ny+2;
  swap->copy0ixa = 1;
  swap->copy0ixb = nx+1;

  // southeast sector; send west

  isector = 2;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nx+1;
  swap->recviy = 0;
  swap->recvnum = nyh+1;
  swap->sendix = 1;
  swap->sendiy = 0;
  swap->sendnum = swap->recvnum;
  swap->recvproc = proceast;
  swap->sendproc = procwest;
  swap->stride = 1;
  swap->copy0ixa = -1;
  swap->copy0ixb = -1;

  // northeast sector; send south

  isector = 3;
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nxh-1;
  swap->recviy = ny+1;
  swap->recvnum = nx-nxh+3;
  swap->sendix = nxh-1;
  swap->sendiy = 1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procnorth;
  swap->sendproc = procsouth;
  swap->stride = ny+2;
  swap->copy0ixa = 1;
  swap->copy0ixb = nx+1;

  // northeast sector; send west

  isector = 3;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nx+1;
  swap->recviy = nyh-1;
  swap->recvnum = ny-nyh+3;
  swap->sendix = 1;
  swap->sendiy = nyh-1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = proceast;
  swap->sendproc = procwest;
  swap->stride = 1;
  swap->copy0ixa = -1;
  swap->copy0ixb = -1;

  maxrecv = 0;
  maxsend = 0;
  for (isector=0;isector<nsector;isector++) {
    for (iswap=0;iswap<nswap;iswap++) {
      swap = &swapinfo[isector][iswap];
      maxrecv=MAX(maxrecv,swap->recvnum);
      maxsend=MAX(maxsend,swap->sendnum);
    }
  }

  recvbuf = (int*) memory->smalloc(maxrecv*sizeof(int),"commgrain:recvbuf");
  sendbuf = (int*) memory->smalloc(maxsend*sizeof(int),"commgrain:sendbuf");
}

/* ----------------------------------------------------------------------
   communication of boundary values for one sector
------------------------------------------------------------------------- */

void CommGrain::communicate(int** lattice, const int isector) {
  int iswap,ix,iy,ixa,ixb,istep,ii,i;
  int *recvpnt, *sendpnt, *recvbufpnt, *sendbufpnt;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;

  // Each swap direction is a little different
  // Need to code each swap explicitly

  // Identify source and destination of data for this swap
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  
  istep = swap->stride;
  ix = swap->recvix;
  iy = swap->recviy;
  recvpnt = &lattice[ix][iy];
  ix = swap->sendix;
  iy = swap->sendiy;
  sendpnt = &lattice[ix][iy];

  // Copy corner to ghost corner

  ixa = swap->copy0ixa;
  ixb = swap->copy0ixb;
  iy = swap->sendiy;
  lattice[ixb][iy] = lattice[ixa][iy];

  // Use MPI to pass data only with other processes

  if (swap->sendproc != universe->me) {

    // Pack the send buffer

    ii = 0;
    for (i=0;i<swap->sendnum;i++) {
      sendbuf[i] = sendpnt[ii];
      ii+=istep;
    }
    recvbufpnt = recvbuf;
    sendbufpnt = sendbuf; 
      
    MPI_Irecv(recvbufpnt,swap->recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbufpnt,swap->sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    for (i=0;i<swap->recvnum;i++) {
      recvpnt[ii] = recvbufpnt[i]; 
      ii+=istep;
    }

  // Use direct access to pass data to self

  } else { 
    ii = 0;
    for (i=0;i<swap->recvnum;i++) {
      recvpnt[ii] = sendpnt[ii];
      ii+=istep;
    }
  }

  // Copy ghost corner to corner

  ixa = swap->copy0ixa;
  ixb = swap->copy0ixb;
  iy = swap->recviy;
  lattice[ixa][iy] = lattice[ixb][iy];

  // Identify source and destination of data for this swap

  iswap = 1;
  swap = &swapinfo[isector][iswap];

  istep = swap->stride;
  ix = swap->recvix;
  iy = swap->recviy;
  recvpnt = &lattice[ix][iy];
  ix = swap->sendix;
  iy = swap->sendiy;
  sendpnt = &lattice[ix][iy];

  // Use MPI to pass data only with other processes

  if (swap->sendproc != universe->me) {

    // Data is contiguous, no buffers needed

    recvbufpnt = recvpnt;
    sendbufpnt = sendpnt; 
      
    MPI_Irecv(recvbufpnt,swap->recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbufpnt,swap->sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
  // Use direct access to pass data to self

  } else { 
    for (i=0;i<swap->recvnum;i+=istep) {
      recvpnt[i] = sendpnt[i];
    }
  }

}

/* ----------------------------------------------------------------------
   deallocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommGrain::free_swap()
{
  memory->destroy_2d_T_array(swapinfo);
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommGrain::allocate_swap(const int idim, const int jdim)
{
  // Allocate a 2D array of SwapInfo objects

  memory->create_2d_T_array(swapinfo,idim,jdim,"commgrain:swapinfo");
  maxsector = idim;
  maxswap = jdim;
}
