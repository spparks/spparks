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
			 const int procsouth_in, const int procnorth_in,
			 const int delghost_in, const int dellocal_in)
{
  nx_local = nx_local_in;
  ny_local = ny_local_in;

  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;

  delghost = delghost_in;
  dellocal = dellocal_in;

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;

  // initialize swap parameters and allocate memory

  SwapInfo* swap;
  int iswap,isector;
  int procarray[2][2];
  SwapInfo swapguide[2][2];
  int ii,jj,kk,irecv,isend;

  // isector = 0,1,2,3 refers to SW,NW,SE,NE quadrants, respectively
  // iswap = 0 refers to the north/south swap
  // iswap = 1 refers to the east/west swap 

  nsector = 4;
  nswap = 2;
  allocate_swap(nsector,nswap);

  // Assign processors to different swap directions as [iswap=0,1][dir=0,1]

  procarray[0][0] = procsouth;
  procarray[0][1] = procnorth;
  procarray[1][0] = procwest;
  procarray[1][1] = proceast;

  // swapguide[0][0] stores values for lower side, not communicating

  swapguide[0][0].recvix = 1-delghost;
  swapguide[0][0].recviy = 1-delghost;
  swapguide[0][0].sendix = 1-delghost;
  swapguide[0][0].sendiy = 1-delghost;
  swapguide[0][0].numx = nx_half-1 + 2*delghost;
  swapguide[0][0].numy = ny_half-1 + 2*delghost;;
  swapguide[0][0].recvproc = -1;
  swapguide[0][0].sendproc = -1;
  swapguide[0][0].copyixa = nx_local-delghost+1;
  swapguide[0][0].copyixb = 1-delghost;

  // swapguide[0][1] stores values for lower side, communicating

  swapguide[0][1].recvix = 1-delghost;
  swapguide[0][1].recviy = 1-delghost;
  swapguide[0][1].sendix = nx_local-delghost+1;
  swapguide[0][1].sendiy = ny_local-delghost+1;
  swapguide[0][1].numx = delghost;
  swapguide[0][1].numy = delghost;
  swapguide[0][1].recvproc = 0;
  swapguide[0][1].sendproc = 1;
  swapguide[0][1].copyixa = -1;
  swapguide[0][1].copyixb = -1;

  // swapguide[1][0] stores values for upper side, not communicating

  swapguide[1][0].recvix = nx_half-delghost;
  swapguide[1][0].recviy = ny_half-delghost;
  swapguide[1][0].sendix = nx_half-delghost;
  swapguide[1][0].sendiy = ny_half-delghost;
  swapguide[1][0].numx = nx_local-nx_half+1+2*delghost;
  swapguide[1][0].numy = ny_local-ny_half+1+2*delghost;
  swapguide[1][0].recvproc = -1;
  swapguide[1][0].sendproc = -1;
  swapguide[1][0].copyixa = 1;
  swapguide[1][0].copyixb = nx_local+1;

  // swapguide[1][1] stores values for upper side, communicating

  swapguide[1][1].recvix = nx_local+1;
  swapguide[1][1].recviy = ny_local+1;
  swapguide[1][1].sendix = 1;
  swapguide[1][1].sendiy = 1;
  swapguide[1][1].numx = delghost;
  swapguide[1][1].numy = delghost;
  swapguide[1][1].recvproc = 1;
  swapguide[1][1].sendproc = 0;
  swapguide[1][1].copyixa = -1;
  swapguide[1][1].copyixb = -1;

  // loop over both sides in each direction, and all swap directions
  // use appropriate guide values in each case

  isector = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      for (iswap = 0; iswap < nswap; iswap++) {
	swap = &swapinfo[isector][iswap];
	
	// identify which direction is communicating	  
	
	if (iswap == 0) {
	  ii = 0;
	  jj = 1;
	  irecv = swapguide[j][jj].recvproc;
	  isend = swapguide[j][jj].sendproc;
	  swap->stride = ny_local+2*delghost;
	} else {
	  ii = 1;
	  jj = 0;
	  irecv = swapguide[i][ii].recvproc;
	  isend = swapguide[i][ii].sendproc;
	  swap->stride = 1;
	}

	swap->recvix = swapguide[i][ii].recvix;
	swap->recviy = swapguide[j][jj].recviy;
	swap->sendix = swapguide[i][ii].sendix;
	swap->sendiy = swapguide[j][jj].sendiy;
	swap->numx = swapguide[i][ii].numx;
	swap->numy = swapguide[j][jj].numy;
	swap->recvproc = procarray[iswap][irecv];
	swap->sendproc = procarray[iswap][isend];
	swap->copyixa = swapguide[i][ii].copyixa;;
	swap->copyixb = swapguide[i][ii].copyixb;
      }
      isector++;
    }

  maxbuf = (MAX(nx_local,ny_local)+2*delghost)*delghost;
  recvbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:recvbuf");
  sendbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:sendbuf");
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
------------------------------------------------------------------------- */

void CommLattice2d::sector(int** lattice, const int isector) {
  int iswap,ii,ix,iy,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;
  int ir,is;

  // Each swap direction is a little different
  // Need to code each swap explicitly

  // First swap is in the y direction

  iswap = 0;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
  int 
    numx = swap->numx,
    numy = swap->numy,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    ixa = swap->copyixa,
    ixb = swap->copyixb;

  recvnum = numx*numy;
  sendnum = recvnum;

  // Copy corner to ghost corner

  ir = ixb;
  is = ixa;
  for (int i=0;i<delghost;i++) {
    iy = sendiy;
    for (int j=0;j<delghost;j++) {
      lattice[ir][iy] = lattice[is][iy]; 
      iy++;
    }
    ir++;
    is++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ix = sendix;
    ii = 0;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	sendbuf[ii] = lattice[ix][iy];
	ii++;
	iy++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ix = recvix;
    ii = 0;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	lattice[ix][iy] = recvbuf[ii]; 
	ii++;
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (int i=0;i<numx;i++) {
      ir = recviy;
      is = sendiy;
      for (int j=0;j<numy;j++) {
	lattice[ix][ir] = lattice[ix][is]; 
	ir++;
	is++;
      }
      ix++;
    }
  }

  // Copy corner to ghost corner

  ir = ixa;
  is = ixb;
  for (int i=0;i<delghost;i++) {
    iy = recviy;
    for (int j=0;j<delghost;j++) {
      lattice[ir][iy] = lattice[is][iy]; 
      iy++;
    }
    ir++;
    is++;
  }


  // Second swap is in the x direction

  iswap = 1;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    ixa = swap->copyixa,
    ixb = swap->copyixb;

  recvnum = numx*numy;
  sendnum = recvnum;

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ix = sendix;
    ii = 0;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	sendbuf[ii] = lattice[ix][iy];
	ii++;
	iy++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ix = recvix;
    ii = 0;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	lattice[ix][iy] = recvbuf[ii];
	ii++;
	iy++;
      }
      ix++;
    }

    // Use direct access to pass data to self

  } else { 

    ir = recvix;
    is = sendix;
//     printf("numx = %d recvix = %d sendix = %d delghost = %d \n ",numx,recvix,sendix,delghost);
//     printf("numy = %d recviy = %d sendiy = %d delghost = %d \n ",numy,recviy,sendiy,delghost);
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	lattice[ir][iy] = lattice[is][iy]; 
	iy++;
      }
      ir++;
      is++;
    }

  }
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
------------------------------------------------------------------------- */

void CommLattice2d::all(int **lattice)
{
  int i,j,m;
  int nsend;
  MPI_Request request;
  MPI_Status status;

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    nsend = nx_local*delghost;
    m = 0;
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= delghost; j++) 
	sendbuf[m++] = lattice[i][j];
    MPI_Irecv(recvbuf,nsend,MPI_INT,procnorth,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= delghost; j++) 
	lattice[i][ny_local+j] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= delghost; j++) 
	sendbuf[m++] = lattice[i][ny_local-delghost+j];
    MPI_Irecv(recvbuf,nsend,MPI_INT,procsouth,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= delghost; j++) 
	lattice[i][j-delghost] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++) 
      for (j = 1; j <= delghost; j++) { 
      lattice[i][ny_local+j] = lattice[i][j];
      lattice[i][j-delghost] = lattice[i][ny_local-delghost+j];
      }
  }
  
  // send west, then east
  // if only 1 proc in x, copy data
  // data is contiguous, so no need for pack/unpack

  if (procwest != me) {
    nsend = (ny_local+2*delghost)*delghost;
    m = 0;
    for (i = 1; i <= delghost; i++) 
      for (j = 1-delghost; j <= ny_local+delghost; j++) 
	sendbuf[m++] = lattice[i][j];
    MPI_Irecv(recvbuf,nsend,MPI_INT,
	      proceast,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procwest,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= delghost; i++) 
      for (j = 1-delghost; j <= ny_local+delghost; j++) 
	lattice[nx_local+i][j] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= delghost; i++) 
      for (j = 1-delghost; j <= ny_local+delghost; j++) 
	sendbuf[m++] = lattice[nx_local-delghost+i][j];
    MPI_Irecv(recvbuf,nsend,MPI_INT,procwest,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,proceast,0,world);
    MPI_Wait(&request,&status);
    m = 0;
    for (i = 1; i <= delghost; i++) 
      for (j = 1-delghost; j <= ny_local+delghost; j++) 
	lattice[i-delghost][j] = recvbuf[m++];


  } else {
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++) {
	lattice[nx_local+i][j] = lattice[i][j];
	lattice[i-delghost][j] = lattice[nx_local-delghost+i][j];
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
