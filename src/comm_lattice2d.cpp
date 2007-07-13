/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice2d.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

CommLattice2d::CommLattice2d(class SPK *spk) : SysPtr(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  swapinfo = NULL;
}

/* ---------------------------------------------------------------------- */

CommLattice2d::~CommLattice2d()
{
  free_swap();
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommLattice2d::init(const int nx_local_in, const int ny_local_in, 
			 const int procwest_in, const int proceast_in, 
			 const int procsouth_in, const int procnorth_in)
{
  nx_local = nx_local_in;
  ny_local = ny_local_in;

  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;

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
  swap->recvnum = nx_half+1;
  swap->sendix = 0;
  swap->sendiy = ny_local;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procsouth;
  swap->sendproc = procnorth;
  swap->stride = ny_local+2;
  swap->copy0ixa = nx_local;
  swap->copy0ixb = 0;

  // southwest sector; send east

  isector = 0;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = 0;
  swap->recvnum = ny_half+1;
  swap->sendix = nx_local;
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
  swap->recviy = ny_local+1;
  swap->recvnum = nx_half+1;
  swap->sendix = 0;
  swap->sendiy = 1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procnorth;
  swap->sendproc = procsouth;
  swap->stride = ny_local+2;
  swap->copy0ixa = nx_local;
  swap->copy0ixb = 0;

  // northwest sector; send east

  isector = 1;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = 0;
  swap->recviy = ny_half-1;
  swap->recvnum = ny_local-ny_half+3;
  swap->sendix = nx_local;
  swap->sendiy = ny_half-1;
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
  swap->recvix = nx_half-1;
  swap->recviy = 0;
  swap->recvnum = nx_local-nx_half+3;
  swap->sendix = nx_half-1;
  swap->sendiy = ny_local;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procsouth;
  swap->sendproc = procnorth;
  swap->stride = ny_local+2;
  swap->copy0ixa = 1;
  swap->copy0ixb = nx_local+1;

  // southeast sector; send west

  isector = 2;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nx_local+1;
  swap->recviy = 0;
  swap->recvnum = ny_half+1;
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
  swap->recvix = nx_half-1;
  swap->recviy = ny_local+1;
  swap->recvnum = nx_local-nx_half+3;
  swap->sendix = nx_half-1;
  swap->sendiy = 1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = procnorth;
  swap->sendproc = procsouth;
  swap->stride = ny_local+2;
  swap->copy0ixa = 1;
  swap->copy0ixb = nx_local+1;

  // northeast sector; send west

  isector = 3;
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  swap->recvix = nx_local+1;
  swap->recviy = ny_half-1;
  swap->recvnum = ny_local-ny_half+3;
  swap->sendix = 1;
  swap->sendiy = ny_half-1;
  swap->sendnum = swap->recvnum;
  swap->recvproc = proceast;
  swap->sendproc = procwest;
  swap->stride = 1;
  swap->copy0ixa = -1;
  swap->copy0ixb = -1;

  maxbuf = MAX(nx_local+2,ny_local+2);
  recvbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:recvbuf");
  sendbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:sendbuf");
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
------------------------------------------------------------------------- */

void CommLattice2d::sector(int** lattice, const int isector) {
  int iswap,ix,iy,ixa,ixb,istep,ii,i;
  int *recvpnt, *sendpnt, *recvbufpnt, *sendbufpnt;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;

  iswap = 0;
  swap = &swapinfo[isector][iswap];
  
  istep = swap->stride;
  ix = swap->recvix;
  iy = swap->recviy;
  recvpnt = &lattice[ix][iy];
  ix = swap->sendix;
  iy = swap->sendiy;
  sendpnt = &lattice[ix][iy];

  // copy corner to ghost corner

  ixa = swap->copy0ixa;
  ixb = swap->copy0ixb;
  iy = swap->sendiy;
  lattice[ixb][iy] = lattice[ixa][iy];

  if (swap->sendproc != me) {
    ii = 0;
    for (i = 0; i < swap->sendnum; i++) {
      sendbuf[i] = sendpnt[ii];
      ii += istep;
    }
    recvbufpnt = recvbuf;
    sendbufpnt = sendbuf; 
      
    MPI_Irecv(recvbufpnt,swap->recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbufpnt,swap->sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
    ii = 0;
    for (i = 0; i < swap->recvnum; i++) {
      recvpnt[ii] = recvbufpnt[i]; 
      ii += istep;
    }

  } else { 
    ii = 0;
    for (i = 0; i < swap->recvnum; i++) {
      recvpnt[ii] = sendpnt[ii];
      ii += istep;
    }
  }

  // copy ghost corner to corner

  ixa = swap->copy0ixa;
  ixb = swap->copy0ixb;
  iy = swap->recviy;
  lattice[ixa][iy] = lattice[ixb][iy];

  iswap = 1;
  swap = &swapinfo[isector][iswap];

  istep = swap->stride;
  ix = swap->recvix;
  iy = swap->recviy;
  recvpnt = &lattice[ix][iy];
  ix = swap->sendix;
  iy = swap->sendiy;
  sendpnt = &lattice[ix][iy];

  if (swap->sendproc != me) {
    recvbufpnt = recvpnt;
    sendbufpnt = sendpnt; 
      
    MPI_Irecv(recvbufpnt,swap->recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbufpnt,swap->sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);

  } else { 
    for (i = 0; i < swap->recvnum; i += istep) {
      recvpnt[i] = sendpnt[i];
    }
  }
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
------------------------------------------------------------------------- */

void CommLattice2d::all(int **lattice)
{
  int i,j,m;
  MPI_Request request;
  MPI_Status status;

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    m = 0;
    for (i = 1; i <= nx_local; i++) sendbuf[m++] = lattice[i][1];
    MPI_Irecv(recvbuf,nx_local,MPI_INT,procnorth,0,world,&request);
    MPI_Send(sendbuf,nx_local,MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= nx_local; i++) lattice[i][ny_local+1] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++) sendbuf[m++] = lattice[i][ny_local];
    MPI_Irecv(recvbuf,nx_local,MPI_INT,procsouth,0,world,&request);
    MPI_Send(sendbuf,nx_local,MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= nx_local; i++) lattice[i][0] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++) {
      lattice[i][ny_local+1] = lattice[i][1];
      lattice[i][0] = lattice[i][ny_local];
    }
  }

  // send west, then east
  // if only 1 proc in x, copy data
  // data is contiguous, so no need for pack/unpack

  if (procwest != me) {
    MPI_Irecv(&lattice[nx_local+1][0],ny_local+2,MPI_INT,
	      proceast,0,world,&request);
    MPI_Send(&lattice[1][0],ny_local+2,MPI_INT,procwest,0,world);
    MPI_Wait(&request,&status);

    MPI_Irecv(&lattice[0][0],ny_local+2,MPI_INT,procwest,0,world,&request);
    MPI_Send(&lattice[nx_local][0],ny_local+2,MPI_INT,proceast,0,world);
    MPI_Wait(&request,&status);

  } else {
    for (j = 0; j <= ny_local+1; j++) {
      lattice[nx_local+1][j] = lattice[1][j];
      lattice[0][j] = lattice[nx_local][j];
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice2d::free_swap()
{
  memory->destroy_2d_T_array(swapinfo);
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice2d::allocate_swap(const int idim, const int jdim)
{
  memory->create_2d_T_array(swapinfo,idim,jdim,"commlattice:swapinfo");
  maxsector = idim;
  maxswap = jdim;
}
