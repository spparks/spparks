/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

CommLattice3d::CommLattice3d(class SPK *spk) : SysPtr(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  swapinfo = NULL;
  reverseinfo = NULL;
  sendbuf = NULL;
  recvbuf = NULL;
}

/* ---------------------------------------------------------------------- */

CommLattice3d::~CommLattice3d()
{
  memory->destroy_2d_T_array(swapinfo);
  memory->destroy_2d_T_array(reverseinfo);
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommLattice3d::init(const int nx_local_in, const int ny_local_in,
			 const int nz_local_in,
			 const int procwest_in, const int proceast_in, 
			 const int procsouth_in, const int procnorth_in,
			 const int procdown_in, const int procup_in,
			 const int delghost_in, const int dellocal_in)
{
  int ii,jj,kk,irecv,isend;

  nx_local = nx_local_in;
  ny_local = ny_local_in;
  nz_local = nz_local_in;

  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;
  procdown = procdown_in;
  procup = procup_in;

  delghost = delghost_in;
  dellocal = dellocal_in;

  // initialize swap parameters and allocate memory

  memory->destroy_2d_T_array(swapinfo);
  memory->destroy_2d_T_array(reverseinfo);
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);

  nsector = 8;
  nswap = 3;
  allocate_swap(nsector,nswap);

  setup_swapinfo();
  setup_reverseinfo();

  maxbuf = (nx_local+2*delghost) * (ny_local*2*delghost) * delghost;
  maxbuf = MAX(maxbuf,(nx_local+2*delghost) * (nz_local*2*delghost) * delghost);
  maxbuf = MAX(maxbuf,(ny_local+2*delghost) * (nz_local*2*delghost) * delghost);
  recvbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:recvbuf");
  sendbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:sendbuf");
}

/* ----------------------------------------------------------------------
   setup up forward communication for each sector
------------------------------------------------------------------------- */

void CommLattice3d::setup_swapinfo()
{
  SwapInfo* swap;
  int iswap,isector;
  int procarray[3][2];
  SwapInfo swapguide[2][2];

  int ii,jj,kk,irecv,isend;

  // initialize swap parameters and allocate memory

  // isector = 0,1,2,3 refer to lower SW,NW,SE,NE quadrants, respectively
  // isector = 4,5,6,7 refer to upper SW,NW,SE,NE quadrants, respectively
  // iswap = 0 refers to the up/down swap
  // iswap = 1 refers to the north/south swap
  // iswap = 2 refers to the east/west swap

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

  // Assign processors to different swap directions as [iswap=0,1,2][dir=0,1]

  procarray[0][0] = procdown;
  procarray[0][1] = procup;
  procarray[1][0] = procsouth;
  procarray[1][1] = procnorth;
  procarray[2][0] = procwest;
  procarray[2][1] = proceast;

  // swapguide[0][0] stores values for lower side, not communicating

  swapguide[0][0].recvix = 1-delghost;
  swapguide[0][0].recviy = 1-delghost;
  swapguide[0][0].recviz = 1-delghost;
  swapguide[0][0].sendix = 1-delghost;
  swapguide[0][0].sendiy = 1-delghost;
  swapguide[0][0].sendiz = 1-delghost;
  swapguide[0][0].numx = nx_half-1 + 2*delghost;;
  swapguide[0][0].numy = ny_half-1 + 2*delghost;;
  swapguide[0][0].numz = nz_half-1 + 2*delghost;;
  swapguide[0][0].recvproc = -1;
  swapguide[0][0].sendproc = -1;
  swapguide[0][0].copyixa = nx_local-delghost+1;
  swapguide[0][0].copyixb = 1-delghost;
  swapguide[0][0].copyiya = ny_local-delghost+1;
  swapguide[0][0].copyiyb = 1-delghost;

  // swapguide[0][1] stores values for lower side, communicating

  swapguide[0][1].recvix = 1-delghost;
  swapguide[0][1].recviy = 1-delghost;
  swapguide[0][1].recviz = 1-delghost;
  swapguide[0][1].sendix = nx_local-delghost+1;
  swapguide[0][1].sendiy = ny_local-delghost+1;
  swapguide[0][1].sendiz = nz_local-delghost+1;
  swapguide[0][1].numx = delghost;
  swapguide[0][1].numy = delghost;
  swapguide[0][1].numz = delghost;
  swapguide[0][1].recvproc = 0;
  swapguide[0][1].sendproc = 1;
  swapguide[0][1].copyixa = -1;
  swapguide[0][1].copyixb = -1;
  swapguide[0][1].copyiya = -1;
  swapguide[0][1].copyiyb = -1;

  // swapguide[1][0] stores values for upper side, not communicating

  swapguide[1][0].recvix = nx_half-delghost;
  swapguide[1][0].recviy = ny_half-delghost;
  swapguide[1][0].recviz = nz_half-delghost;
  swapguide[1][0].sendix = nx_half-delghost;
  swapguide[1][0].sendiy = ny_half-delghost;
  swapguide[1][0].sendiz = nz_half-delghost;
  swapguide[1][0].numx = nx_local-nx_half+1+2*delghost;
  swapguide[1][0].numy = ny_local-ny_half+1+2*delghost;
  swapguide[1][0].numz = nz_local-nz_half+1+2*delghost;
  swapguide[1][0].recvproc = -1;
  swapguide[1][0].sendproc = -1;
  swapguide[1][0].copyixa = 1;
  swapguide[1][0].copyixb = nx_local+1;
  swapguide[1][0].copyiya = 1;
  swapguide[1][0].copyiyb = ny_local+1;

  // swapguide[1][1] stores values for upper side, communicating

  swapguide[1][1].recvix = nx_local+1;
  swapguide[1][1].recviy = ny_local+1;
  swapguide[1][1].recviz = nz_local+1;
  swapguide[1][1].sendix = 1;
  swapguide[1][1].sendiy = 1;
  swapguide[1][1].sendiz = 1;
  swapguide[1][1].numx = delghost;
  swapguide[1][1].numy = delghost;
  swapguide[1][1].numz = delghost;
  swapguide[1][1].recvproc = 1;
  swapguide[1][1].sendproc = 0;
  swapguide[1][1].copyixa = -1;
  swapguide[1][1].copyixb = -1;
  swapguide[1][1].copyiya = -1;
  swapguide[1][1].copyiyb = -1;

  // loop over both sides in each direction, and all swap directions
  // use appropriate guide values in each case

  isector = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	for (iswap = 0; iswap < nswap; iswap++) {
	  swap = &swapinfo[isector][iswap];

	  // identify which direction is communicating	  

	  if (iswap == 0) {
	    ii = 0;
	    jj = 0;
	    kk = 1;
	    irecv = swapguide[k][kk].recvproc;
	    isend = swapguide[k][kk].sendproc;
	  }  else if (iswap == 1) {
	    ii = 0;
	    jj = 1;
	    kk = 0;
	    irecv = swapguide[j][jj].recvproc;
	    isend = swapguide[j][jj].sendproc;
	  } else {
	    ii = 1;
	    jj = 0;
	    kk = 0;
	    irecv = swapguide[i][ii].recvproc;
	    isend = swapguide[i][ii].sendproc;
	  }

	  swap->recvix = swapguide[i][ii].recvix;
	  swap->recviy = swapguide[j][jj].recviy;
	  swap->recviz = swapguide[k][kk].recviz;
	  swap->sendix = swapguide[i][ii].sendix;
	  swap->sendiy = swapguide[j][jj].sendiy;
	  swap->sendiz = swapguide[k][kk].sendiz;
	  swap->numx = swapguide[i][ii].numx;
	  swap->numy = swapguide[j][jj].numy;
	  swap->numz = swapguide[k][kk].numz;
	  swap->recvproc = procarray[iswap][irecv];
	  swap->sendproc = procarray[iswap][isend];
	  swap->copyixa = swapguide[i][ii].copyixa;;
	  swap->copyixb = swapguide[i][ii].copyixb;
	  swap->copyiya = swapguide[j][jj].copyiya;
	  swap->copyiyb = swapguide[j][jj].copyiyb;
	}
	isector++;
      }

}

/* ----------------------------------------------------------------------
   setup up reverse communication for each sector
------------------------------------------------------------------------- */

void CommLattice3d::setup_reverseinfo()
{
  int ii,jj,kk,irecv,isend;

  // initialize swap parameters and allocate memory

  SwapInfo* swap;
  int iswap,isector;
  int procarray[3][2];
  SwapInfo swapguide[2][2];

  // isector = 0,1,2,3 refer to lower SW,NW,SE,NE quadrants, respectively
  // isector = 4,5,6,7 refer to upper SW,NW,SE,NE quadrants, respectively
  // iswap = 0 refers to the up/down swap
  // iswap = 1 refers to the north/south swap
  // iswap = 2 refers to the east/west swap

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

  // Assign processors to different swap directions as [iswap=0,1,2][dir=0,1]

  procarray[0][0] = procdown;
  procarray[0][1] = procup;
  procarray[1][0] = procsouth;
  procarray[1][1] = procnorth;
  procarray[2][0] = procwest;
  procarray[2][1] = proceast;

  // swapguide[0][0] stores values for lower side, not communicating

  swapguide[0][0].recvix = 1-dellocal;
  swapguide[0][0].recviy = 1-dellocal;
  swapguide[0][0].recviz = 1-dellocal;
  swapguide[0][0].sendix = 1-dellocal;
  swapguide[0][0].sendiy = 1-dellocal;
  swapguide[0][0].sendiz = 1-dellocal;
  swapguide[0][0].numx = nx_half-1 + 2*dellocal;;
  swapguide[0][0].numy = ny_half-1 + 2*dellocal;;
  swapguide[0][0].numz = nz_half-1 + 2*dellocal;;
  swapguide[0][0].recvproc = -1;
  swapguide[0][0].sendproc = -1;
  swapguide[0][0].copyixa = nx_local-dellocal+1;
  swapguide[0][0].copyixb = 1-dellocal;
  swapguide[0][0].copyiya = ny_local-dellocal+1;
  swapguide[0][0].copyiyb = 1-dellocal;

  // swapguide[0][1] stores values for lower side, communicating

  swapguide[0][1].recvix = 1-dellocal;
  swapguide[0][1].recviy = 1-dellocal;
  swapguide[0][1].recviz = 1-dellocal;
  swapguide[0][1].sendix = nx_local-dellocal+1;
  swapguide[0][1].sendiy = ny_local-dellocal+1;
  swapguide[0][1].sendiz = nz_local-dellocal+1;
  swapguide[0][1].numx = dellocal;
  swapguide[0][1].numy = dellocal;
  swapguide[0][1].numz = dellocal;
  swapguide[0][1].recvproc = 0;
  swapguide[0][1].sendproc = 1;
  swapguide[0][1].copyixa = -1;
  swapguide[0][1].copyixb = -1;
  swapguide[0][1].copyiya = -1;
  swapguide[0][1].copyiyb = -1;

  // swapguide[1][0] stores values for upper side, not communicating

  swapguide[1][0].recvix = nx_half-dellocal;
  swapguide[1][0].recviy = ny_half-dellocal;
  swapguide[1][0].recviz = nz_half-dellocal;
  swapguide[1][0].sendix = nx_half-dellocal;
  swapguide[1][0].sendiy = ny_half-dellocal;
  swapguide[1][0].sendiz = nz_half-dellocal;
  swapguide[1][0].numx = nx_local-nx_half+1+2*dellocal;
  swapguide[1][0].numy = ny_local-ny_half+1+2*dellocal;
  swapguide[1][0].numz = nz_local-nz_half+1+2*dellocal;
  swapguide[1][0].recvproc = -1;
  swapguide[1][0].sendproc = -1;
  swapguide[1][0].copyixa = 1;
  swapguide[1][0].copyixb = nx_local+1;
  swapguide[1][0].copyiya = 1;
  swapguide[1][0].copyiyb = ny_local+1;

  // swapguide[1][1] stores values for upper side, communicating

  swapguide[1][1].recvix = nx_local+1;
  swapguide[1][1].recviy = ny_local+1;
  swapguide[1][1].recviz = nz_local+1;
  swapguide[1][1].sendix = 1;
  swapguide[1][1].sendiy = 1;
  swapguide[1][1].sendiz = 1;
  swapguide[1][1].numx = dellocal;
  swapguide[1][1].numy = dellocal;
  swapguide[1][1].numz = dellocal;
  swapguide[1][1].recvproc = 1;
  swapguide[1][1].sendproc = 0;
  swapguide[1][1].copyixa = -1;
  swapguide[1][1].copyixb = -1;
  swapguide[1][1].copyiya = -1;
  swapguide[1][1].copyiyb = -1;

  // loop over both sides in each direction, and all swap directions
  // use appropriate guide values in each case

  isector = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	for (iswap = 0; iswap < nswap; iswap++) {
	  // Same order as setup_swapinfo(); order reversed in reverse_sector();
	  swap = &reverseinfo[isector][iswap];

	  // identify which direction is communicating	  

	  if (iswap == 0) {
	    ii = 0;
	    jj = 0;
	    kk = 1;
	    irecv = swapguide[k][kk].recvproc;
	    isend = swapguide[k][kk].sendproc;
	  }  else if (iswap == 1) {
	    ii = 0;
	    jj = 1;
	    kk = 0;
	    irecv = swapguide[j][jj].recvproc;
	    isend = swapguide[j][jj].sendproc;
	  } else {
	    ii = 1;
	    jj = 0;
	    kk = 0;
	    irecv = swapguide[i][ii].recvproc;
	    isend = swapguide[i][ii].sendproc;
	  }

	  // Exactly opposite to assignments in setup_swapinfo()
	  swap->recvix = swapguide[i][ii].sendix;
	  swap->recviy = swapguide[j][jj].sendiy;
	  swap->recviz = swapguide[k][kk].sendiz;
	  swap->sendix = swapguide[i][ii].recvix;
	  swap->sendiy = swapguide[j][jj].recviy;
	  swap->sendiz = swapguide[k][kk].recviz;
	  swap->numx = swapguide[i][ii].numx;
	  swap->numy = swapguide[j][jj].numy;
	  swap->numz = swapguide[k][kk].numz;
	  swap->recvproc = procarray[iswap][isend];
	  swap->sendproc = procarray[iswap][irecv];
	  // These values are not reversed; it's complicated
	  swap->copyixa = swapguide[i][ii].copyixa;
	  swap->copyixb = swapguide[i][ii].copyixb;
	  swap->copyiya = swapguide[j][jj].copyiya;
	  swap->copyiyb = swapguide[j][jj].copyiyb;
	}
	isector++;
      }

}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   choose one of two different variants
------------------------------------------------------------------------- */

void CommLattice3d::sector(int*** lattice, const int isector) {
  if (delghost == 1) sector_onelayer(lattice,isector);
  else sector_multilayer(lattice,isector);
}

// // This was a version used for testing
// void CommLattice3d::sector(int*** lattice, const int isector) {
//   if (delghost == 1) {
//     // This is just a test
//     sector_onelayer(lattice,isector);
//   }
//   else {
//     // This is just a test
//     sector_multilayer_destroy(lattice,isector);
//     reverse_sector_multilayer(lattice,isector);
//     sector_multilayer(lattice,isector);
//   }
// }

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   delghost = 1
------------------------------------------------------------------------- */

void CommLattice3d::sector_onelayer(int*** lattice, const int isector) {
  int i,j,iswap,ii,ix,iy,iz,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;
  int ir,is,jr,js,kr,ks;

  // Each swap direction is a little different
  // Need to code each swap explicitly

  // First swap is in the z direction

  iswap = 0;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
  int 
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numy;
  sendnum = recvnum;

  // Copy y-edge to ghost edge
  iy = sendiy;
  for (i=0;i<numy;i++) {
    lattice[ixb][iy][sendiz] = lattice[ixa][iy][sendiz];
    iy++;
  }

  // Copy x-edge to ghost edge

  ix = sendix;
  for (i=0;i<numx;i++) {
    lattice[ix][iyb][sendiz] = lattice[ix][iya][sendiz];
    ix++;
  }

  // Copy corner to ghost corner; must do this last to avoid overwrite

  lattice[ixb][iyb][sendiz] = lattice[ixa][iya][sendiz];

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (i=0;i<numx;i++) {
      iy = sendiy;
      for (j=0;j<numy;j++) {
	sendbuf[ii] = lattice[ix][iy][sendiz];
	iy++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    ix = recvix;
    for (i=0;i<numx;i++) {
      iy = recviy;
      for (j=0;j<numy;j++) {
	lattice[ix][iy][recviz] = recvbuf[ii]; 
	iy++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (i=0;i<numx;i++) {
      iy = recviy;
      for (j=0;j<numy;j++) {
	lattice[ix][iy][recviz] = lattice[ix][iy][sendiz]; 
	iy++;
      }
      ix++;
    }
  }

  // Copy ghost y-edge to edge

  iy = recviy;
  for (i=0;i<numy;i++) {
    lattice[ixa][iy][recviz] = lattice[ixb][iy][recviz];
    iy++;
  }

  // Copy ghost x-edge to edge

  ix = recvix;
  for (i=0;i<numx;i++) {
    lattice[ix][iya][recviz] = lattice[ix][iyb][recviz];
    ix++;
  }

  // Copy ghost corner to corner; must do this last to avoid overwrite

  lattice[ixa][iya][recviz] = lattice[ixb][iyb][recviz];

  // Second swap is in the y direction

  iswap = 1;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numz;
  sendnum = recvnum;
  
  // Copy edge to ghost edge

  iz = sendiz;
  for (i=0;i<numz;i++) {
    lattice[ixb][sendiy][iz] = lattice[ixa][sendiy][iz];
    iz++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (i=0;i<numx;i++) {
      iz = sendiz;
      for (j=0;j<numz;j++) {
	sendbuf[ii] = lattice[ix][sendiy][iz];
	iz++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    ix = recvix;
    for (i=0;i<numx;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[ix][recviy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (i=0;i<numx;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[ix][recviy][iz] = lattice[ix][sendiy][iz]; 
	iz++;
      }
      ix++;
    }
  }

  // Copy ghost edge to edge

  iz = recviz;
  for (i=0;i<numz;i++) {
    lattice[ixa][recviy][iz] = lattice[ixb][recviy][iz];
    iz++;
  }

  // Third swap is in the x direction

  iswap = 2;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numy*numz;
  sendnum = recvnum;

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    iy = sendiy;
    for (i=0;i<numy;i++) {
      iz = sendiz;
      for (j=0;j<numz;j++) {
	sendbuf[ii] = lattice[sendix][iy][iz];
	iz++;
	ii++;
      }
      iy++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    iy = recviy;
    for (i=0;i<numy;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[recvix][iy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      iy++;
    }

  // Use direct access to pass data to self

  } else { 

    iy = recviy;
    for (i=0;i<numy;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[recvix][iy][iz] = lattice[sendix][iy][iz]; 
	iz++;
      }
      iy++;
    }
  }
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   delghost > 1
------------------------------------------------------------------------- */

void CommLattice3d::sector_multilayer(int*** lattice, const int isector) {
  int iswap,ii,ix,iy,iz,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;
  int ir,is,jr,js,kr,ks;

  // Each swap direction is a little different
  // Need to code each swap explicitly

  // First swap is in the z direction

  iswap = 0;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
  int 
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numy*delghost;
  sendnum = recvnum;

  // Copy y-edge to ghost edge
  ir = ixb;
  is = ixa;
  for (int i=0;i<delghost;i++) {
    iy = sendiy;
    for (int j=0;j<numy;j++) {
      iz = sendiz;
      for (int k=0;k<delghost;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Copy x-edge to ghost edge

  ix = sendix;
  for (int i=0;i<numx;i++) {
    jr = iyb;
    js = iya;
    for (int j=0;j<delghost;j++) {
      iz = sendiz;
      for (int k=0;k<delghost;k++) {
	lattice[ix][jr][iz] = lattice[ix][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ix++;
  }

  // Copy corner to ghost corner; must do this last to avoid overwrite

  ir = ixb;
  is = ixa;
  for (int i=0;i<delghost;i++) {
    jr = iyb;
    js = iya;
    for (int j=0;j<delghost;j++) {
      iz = sendiz;
      for (int k=0;k<delghost;k++) {
	lattice[ir][jr][iz] = lattice[is][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ir++;
    is++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	iz = sendiz;
	for (int k=0;k<delghost;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<delghost;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	kr = recviz;
	ks = sendiz;
	for (int k=0;k<delghost;k++) {
	  lattice[ix][iy][kr] = lattice[ix][iy][ks]; 
	  kr++;
	  ks++;
	}
	iy++;
      }
      ix++;
    }

  }

  // Copy ghost y-edge to edge

  ir = ixa;
  is = ixb;
  for (int i=0;i<delghost;i++) {
    iy = recviy;
    for (int j=0;j<numy;j++) {
      iz = recviz;
      for (int k=0;k<delghost;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Copy ghost x-edge to edge

  ix = recvix;
  for (int i=0;i<numx;i++) {
    jr = iya;
    js = iyb;
    for (int j=0;j<delghost;j++) {
      iz = recviz;
      for (int k=0;k<delghost;k++) {
	lattice[ix][jr][iz] = lattice[ix][js][iz];
	iz++;
      }
      jr++;
      js++;
    }
    ix++;
  }

  // Copy ghost corner to corner; must do this last to avoid overwrite

  ir = ixa;
  is = ixb;
  for (int i=0;i<delghost;i++) {
    jr = iya;
    js = iyb;
    for (int j=0;j<delghost;j++) {
      iz = recviz;
      for (int k=0;k<delghost;k++) {
	lattice[ir][jr][iz] = lattice[is][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ir++;
    is++;
  }

  // Second swap is in the y direction

  iswap = 1;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numz*delghost;
  sendnum = recvnum;
  
  // Copy z-edge to ghost edge

  ir = ixb;
  is = ixa;
  for (int i=0;i<delghost;i++) {
    iy = sendiy;
    for (int j=0;j<delghost;j++) {
      iz = sendiz;
      for (int k=0;k<numz;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<delghost;j++) {
	iz = sendiz;
	for (int k=0;k<numz;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<delghost;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (int i=0;i<numx;i++) {
      jr = recviy;
      js = sendiy;
      for (int j=0;j<delghost;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][jr][iz] = lattice[ix][js][iz]; 
	  iz++;
	}
	jr++;
	js++;
      }
      ix++;
    }

  }

  // Copy ghost edge to z-edge

  ir = ixa;
  is = ixb;
  for (int i=0;i<delghost;i++) {
    iy = recviy;
    for (int j=0;j<delghost;j++) {
      iz = recviz;
      for (int k=0;k<numz;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Third swap is in the x direction

  iswap = 2;
  swap = &swapinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numy*numz*delghost;
  sendnum = recvnum;

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<delghost;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	iz = sendiz;
	for (int k=0;k<numz;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<delghost;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ir = recvix;
    is = sendix;
    for (int i=0;i<delghost;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	  iz++;
	}
      iy++;
      }
      ir++;
      is++;
    }

  }
}


/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   choose one of two different variants
------------------------------------------------------------------------- */

void CommLattice3d::reverse_sector(int*** lattice, const int isector) {
  if (dellocal == 0) return;
  else if (dellocal == 1) reverse_sector_onelayer(lattice,isector);
  else reverse_sector_multilayer(lattice,isector);
}

/* ----------------------------------------------------------------------
   send back ghost values for one sector (quadrant)
   dellocal = 1
------------------------------------------------------------------------- */

void CommLattice3d::reverse_sector_onelayer(int*** lattice, const int isector) {
  int i,j,iswap,ii,ix,iy,iz,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;
  int ir,is,jr,js,kr,ks;

  // Each swap direction is a little different
  // Need to code each swap explicitly
  // These swaps are exactly the opposite of those in sector()

  // Third swap is in the x direction

  iswap = 2;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
  int
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numy*numz;
  sendnum = recvnum;

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    iy = sendiy;
    for (i=0;i<numy;i++) {
      iz = sendiz;
      for (j=0;j<numz;j++) {
	sendbuf[ii] = lattice[sendix][iy][iz];
	iz++;
	ii++;
      }
      iy++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    iy = recviy;
    for (i=0;i<numy;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[recvix][iy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      iy++;
    }

  // Use direct access to pass data to self

  } else { 

    iy = recviy;
    for (i=0;i<numy;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[recvix][iy][iz] = lattice[sendix][iy][iz]; 
	iz++;
      }
      iy++;
    }
  }

  // Second swap is in the y direction

  iswap = 1;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numz;
  sendnum = recvnum;
  
  // Copy edge to ghost edge

  iz = sendiz;
  for (i=0;i<numz;i++) {
    lattice[ixb][sendiy][iz] = lattice[ixa][sendiy][iz];
    iz++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (i=0;i<numx;i++) {
      iz = sendiz;
      for (j=0;j<numz;j++) {
	sendbuf[ii] = lattice[ix][sendiy][iz];
	iz++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    ix = recvix;
    for (i=0;i<numx;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[ix][recviy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (i=0;i<numx;i++) {
      iz = recviz;
      for (j=0;j<numz;j++) {
	lattice[ix][recviy][iz] = lattice[ix][sendiy][iz]; 
	iz++;
      }
      ix++;
    }
  }

  // Copy ghost edge to edge

  iz = recviz;
  for (i=0;i<numz;i++) {
    lattice[ixa][recviy][iz] = lattice[ixb][recviy][iz];
    iz++;
  }

  // Third swap is in the z direction

  iswap = 0;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numy;
  sendnum = recvnum;

  // Copy y-edge to ghost edge
  iy = sendiy;
  for (i=0;i<numy;i++) {
    lattice[ixb][iy][sendiz] = lattice[ixa][iy][sendiz];
    iy++;
  }

  // Copy x-edge to ghost edge

  ix = sendix;
  for (i=0;i<numx;i++) {
    lattice[ix][iyb][sendiz] = lattice[ix][iya][sendiz];
    ix++;
  }

  // Copy corner to ghost corner; must do this last to avoid overwrite

  lattice[ixb][iyb][sendiz] = lattice[ixa][iya][sendiz];

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (i=0;i<numx;i++) {
      iy = sendiy;
      for (j=0;j<numy;j++) {
	sendbuf[ii] = lattice[ix][iy][sendiz];
	iy++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    ix = recvix;
    for (i=0;i<numx;i++) {
      iy = recviy;
      for (j=0;j<numy;j++) {
	lattice[ix][iy][recviz] = recvbuf[ii]; 
	iy++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (i=0;i<numx;i++) {
      iy = recviy;
      for (j=0;j<numy;j++) {
	lattice[ix][iy][recviz] = lattice[ix][iy][sendiz]; 
	iy++;
      }
      ix++;
    }
  }

  // Copy ghost y-edge to edge

  iy = recviy;
  for (i=0;i<numy;i++) {
    lattice[ixa][iy][recviz] = lattice[ixb][iy][recviz];
    iy++;
  }

  // Copy ghost x-edge to edge

  ix = recvix;
  for (i=0;i<numx;i++) {
    lattice[ix][iya][recviz] = lattice[ix][iyb][recviz];
    ix++;
  }

  // Copy ghost corner to corner; must do this last to avoid overwrite

  lattice[ixa][iya][recviz] = lattice[ixb][iyb][recviz];

}

/* ----------------------------------------------------------------------
   send back ghost values for one sector (quadrant)
   dellocal > 1
------------------------------------------------------------------------- */

void CommLattice3d::reverse_sector_multilayer(int*** lattice, const int isector) {
  int iswap,ii,ix,iy,iz,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;
  int ir,is,jr,js,kr,ks;

  // Each swap direction is a little different
  // Need to code each swap explicitly
  // These swaps are exactly the opposite of those in sector()

  // First swap is in the x direction

  iswap = 2;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
  int
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numy*numz*dellocal;
  sendnum = recvnum;

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<dellocal;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	iz = sendiz;
	for (int k=0;k<numz;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<dellocal;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ir = recvix;
    is = sendix;
    for (int i=0;i<dellocal;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	  iz++;
	}
      iy++;
      }
      ir++;
      is++;
    }

  }

  // Second swap is in the y direction

  iswap = 1;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numz*dellocal;
  sendnum = recvnum;
  
  // Copy z-edge to ghost edge

  ir = ixb;
  is = ixa;
  for (int i=0;i<dellocal;i++) {
    iy = sendiy;
    for (int j=0;j<dellocal;j++) {
      iz = sendiz;
      for (int k=0;k<numz;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<dellocal;j++) {
	iz = sendiz;
	for (int k=0;k<numz;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<dellocal;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (int i=0;i<numx;i++) {
      jr = recviy;
      js = sendiy;
      for (int j=0;j<dellocal;j++) {
	iz = recviz;
	for (int k=0;k<numz;k++) {
	  lattice[ix][jr][iz] = lattice[ix][js][iz]; 
	  iz++;
	}
	jr++;
	js++;
      }
      ix++;
    }

  }

  // Copy ghost edge to z-edge

  ir = ixa;
  is = ixb;
  for (int i=0;i<dellocal;i++) {
    iy = recviy;
    for (int j=0;j<dellocal;j++) {
      iz = recviz;
      for (int k=0;k<numz;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // First swap is in the z direction

  iswap = 0;
  swap = &reverseinfo[isector][iswap];

  // local copies of SwapInfo data
    numx = swap->numx,
    numy = swap->numy,
    numz = swap->numz,
    sendproc = swap->sendproc,
    recvproc = swap->recvproc,
    sendix = swap->sendix,
    recvix = swap->recvix,
    sendiy = swap->sendiy,
    recviy = swap->recviy,
    sendiz = swap->sendiz,
    recviz = swap->recviz,
    ixa = swap->copyixa,
    ixb = swap->copyixb,
    iya = swap->copyiya,
    iyb = swap->copyiyb;

  recvnum = numx*numy*dellocal;
  sendnum = recvnum;

  // Copy y-edge to ghost edge
  ir = ixb;
  is = ixa;
  for (int i=0;i<dellocal;i++) {
    iy = sendiy;
    for (int j=0;j<numy;j++) {
      iz = sendiz;
      for (int k=0;k<dellocal;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Copy x-edge to ghost edge

  ix = sendix;
  for (int i=0;i<numx;i++) {
    jr = iyb;
    js = iya;
    for (int j=0;j<dellocal;j++) {
      iz = sendiz;
      for (int k=0;k<dellocal;k++) {
	lattice[ix][jr][iz] = lattice[ix][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ix++;
  }

  // Copy corner to ghost corner; must do this last to avoid overwrite

  ir = ixb;
  is = ixa;
  for (int i=0;i<dellocal;i++) {
    jr = iyb;
    js = iya;
    for (int j=0;j<dellocal;j++) {
      iz = sendiz;
      for (int k=0;k<dellocal;k++) {
	lattice[ir][jr][iz] = lattice[is][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ir++;
    is++;
  }

  // Use MPI to pass data only with other processes

  if (sendproc != me) {

    // Pack the send buffer

    ii = 0;
    ix = sendix;
    for (int i=0;i<numx;i++) {
      iy = sendiy;
      for (int j=0;j<numy;j++) {
	iz = sendiz;
	for (int k=0;k<dellocal;k++) {
	  sendbuf[ii] = lattice[ix][iy][iz];
	  iz++;
	  ii++;
	}
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

    ii = 0;
    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	iz = recviz;
	for (int k=0;k<dellocal;k++) {
	  lattice[ix][iy][iz] = recvbuf[ii]; 
	  iz++;
	  ii++;
	}
	iy++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = recvix;
    for (int i=0;i<numx;i++) {
      iy = recviy;
      for (int j=0;j<numy;j++) {
	kr = recviz;
	ks = sendiz;
	for (int k=0;k<dellocal;k++) {
	  lattice[ix][iy][kr] = lattice[ix][iy][ks]; 
	  kr++;
	  ks++;
	}
	iy++;
      }
      ix++;
    }

  }

  // Copy ghost y-edge to edge

  ir = ixa;
  is = ixb;
  for (int i=0;i<dellocal;i++) {
    iy = recviy;
    for (int j=0;j<numy;j++) {
      iz = recviz;
      for (int k=0;k<dellocal;k++) {
	lattice[ir][iy][iz] = lattice[is][iy][iz]; 
	iz++;
      }
      iy++;
    }
    ir++;
    is++;
  }

  // Copy ghost x-edge to edge

  ix = recvix;
  for (int i=0;i<numx;i++) {
    jr = iya;
    js = iyb;
    for (int j=0;j<dellocal;j++) {
      iz = recviz;
      for (int k=0;k<dellocal;k++) {
	lattice[ix][jr][iz] = lattice[ix][js][iz];
	iz++;
      }
      jr++;
      js++;
    }
    ix++;
  }

  // Copy ghost corner to corner; must do this last to avoid overwrite

  ir = ixa;
  is = ixb;
  for (int i=0;i<dellocal;i++) {
    jr = iya;
    js = iyb;
    for (int j=0;j<dellocal;j++) {
      iz = recviz;
      for (int k=0;k<dellocal;k++) {
	lattice[ir][jr][iz] = lattice[is][js][iz]; 
	iz++;
      }
      jr++;
      js++;
    }
    ir++;
    is++;
  }

}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
   choose one of two different variants
------------------------------------------------------------------------- */

void CommLattice3d::all(int*** lattice) {
  if (delghost == 1) all_onelayer(lattice);
  else all_multilayer(lattice);
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
   delghost = 1
------------------------------------------------------------------------- */

void CommLattice3d::all_onelayer(int ***lattice)
{
  int i,j,k,m;
  MPI_Request request;
  MPI_Status status;

  // send down, then up
  // if only 1 proc in z, copy data

  if (procdown != me) {
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	sendbuf[m++] = lattice[i][j][1];

    MPI_Irecv(recvbuf,nx_local*ny_local,MPI_INT,procup,0,world,&request);
    MPI_Send(sendbuf,nx_local*ny_local,MPI_INT,procdown,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	lattice[i][j][nz_local+1] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	sendbuf[m++] = lattice[i][j][nz_local];

    MPI_Irecv(recvbuf,nx_local*ny_local,MPI_INT,procdown,0,world,&request);
    MPI_Send(sendbuf,nx_local*ny_local,MPI_INT,procup,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	lattice[i][j][0] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++) {
	lattice[i][j][nz_local+1] = lattice[i][j][1];
	lattice[i][j][0] = lattice[i][j][nz_local];
      }
  }

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	sendbuf[m++] = lattice[i][1][k];

    MPI_Irecv(recvbuf,nx_local*(nz_local+2),MPI_INT,procnorth,0,world,
	      &request);
    MPI_Send(sendbuf,nx_local*(nz_local+2),MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	lattice[i][ny_local+1][k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	sendbuf[m++] = lattice[i][ny_local][k];

    MPI_Irecv(recvbuf,nx_local*(nz_local+2),MPI_INT,procsouth,0,world,
	      &request);
    MPI_Send(sendbuf,nx_local*(nz_local+2),MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	lattice[i][0][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++) {
	lattice[i][ny_local+1][k] = lattice[i][1][k];
	lattice[i][0][k] = lattice[i][ny_local][k];
      }
  }

  // send west, then east
  // if only 1 proc in x, copy data
  // data is contiguous, so no need for pack/unpack

  if (procwest != me) {
    MPI_Irecv(&lattice[nx_local+1][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	      proceast,0,world,&request);
    MPI_Send(&lattice[1][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	     procwest,0,world);
    MPI_Wait(&request,&status);

    MPI_Irecv(&lattice[0][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	      procwest,0,world,&request);
    MPI_Send(&lattice[nx_local][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	     proceast,0,world);
    MPI_Wait(&request,&status);

  } else {
    for (j = 0; j <= ny_local+1; j++)
      for (k = 0; k <= nz_local+1; k++) {
	lattice[nx_local+1][j][k] = lattice[1][j][k];
	lattice[0][j][k] = lattice[nx_local][j][k];
      }
  }
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
   delghost > 1
------------------------------------------------------------------------- */

void CommLattice3d::all_multilayer(int ***lattice)
{
  int i,j,k,m,nsend;
  MPI_Request request;
  MPI_Status status;

  // send down, then up
  // if only 1 proc in z, copy data

  if (procdown != me) {
    nsend = nx_local*ny_local*delghost;
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= delghost; k++)
	sendbuf[m++] = lattice[i][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procup,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procdown,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= delghost; k++)
	  lattice[i][j][nz_local+k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= delghost; k++)
	  sendbuf[m++] = lattice[i][j][nz_local-delghost+k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procdown,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procup,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= delghost; k++)
	  lattice[i][j][k-delghost] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= delghost; k++) {
	lattice[i][j][nz_local+k] = lattice[i][j][k];
	lattice[i][j][k-delghost] = lattice[i][j][nz_local-delghost+k];
	}
  }

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    nsend = nx_local*(nz_local+2*delghost)*delghost;
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  sendbuf[m++] = lattice[i][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procnorth,0,world,
	      &request);
    MPI_Send(sendbuf,nsend,MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  lattice[i][ny_local+j][k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  sendbuf[m++] = lattice[i][ny_local-delghost+j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procsouth,0,world,
	      &request);
    MPI_Send(sendbuf,nsend,MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  lattice[i][j-delghost][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++) {
	  lattice[i][ny_local+j][k] = lattice[i][j][k];
	  lattice[i][j-delghost][k] = lattice[i][ny_local-delghost+j][k];
	}
  }

  // send west, then east
  // if only 1 proc in x, copy data

  if (procwest != me) {
    nsend = (ny_local+2*delghost)*(nz_local+2*delghost)*delghost;
    m = 0;
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  sendbuf[m++] = lattice[i][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,proceast,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procwest,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  lattice[nx_local+i][j][k] = recvbuf[m++];


    m = 0;
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  sendbuf[m++] = lattice[nx_local-delghost+i][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procwest,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,proceast,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++)
	  lattice[i-delghost][j][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= delghost; i++)
      for (j = 1-delghost; j <= ny_local+delghost; j++)
	for (k = 1-delghost; k <= nz_local+delghost; k++) {
	  lattice[nx_local+i][j][k] = lattice[i][j][k];
	  lattice[i-delghost][j][k] = lattice[nx_local-delghost+i][j][k];
	}
  }
}

/* ----------------------------------------------------------------------
   send back ghost values for entire sub-domain owned by this proc
   choose one of two different variants
------------------------------------------------------------------------- */

void CommLattice3d::reverse_all(int*** lattice) {
  if (dellocal == 0) return;
  else if (dellocal == 1) reverse_all_onelayer(lattice);
  else reverse_all_multilayer(lattice);
}

/* ----------------------------------------------------------------------
   send back ghost values for entire sub-domain owned by this proc
   dellocal = 1
------------------------------------------------------------------------- */

void CommLattice3d::reverse_all_onelayer(int ***lattice)
{
  int i,j,k,m;
  MPI_Request request;
  MPI_Status status;

  // send west, then east
  // if only 1 proc in x, copy data
  // data is contiguous, so no need for pack/unpack

  if (procwest != me) {
    MPI_Irecv(&lattice[1][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	      proceast,0,world,&request);
    MPI_Send(&lattice[nx_local+1][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	     procwest,0,world);
    MPI_Wait(&request,&status);

    MPI_Irecv(&lattice[nx_local][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	      procwest,0,world,&request);
    MPI_Send(&lattice[0][0][0],(ny_local+2)*(nz_local+2),MPI_INT,
	     proceast,0,world);
    MPI_Wait(&request,&status);

  } else {
    for (j = 0; j <= ny_local+1; j++)
      for (k = 0; k <= nz_local+1; k++) {
	lattice[1][j][k] = lattice[nx_local+1][j][k];
	lattice[nx_local][j][k] = lattice[0][j][k];
      }
  }

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	sendbuf[m++] = lattice[i][ny_local+1][k];

    MPI_Irecv(recvbuf,nx_local*(nz_local+2),MPI_INT,procnorth,0,world,
	      &request);
    MPI_Send(sendbuf,nx_local*(nz_local+2),MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	lattice[i][1][k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	sendbuf[m++] = lattice[i][0][k];

    MPI_Irecv(recvbuf,nx_local*(nz_local+2),MPI_INT,procsouth,0,world,
	      &request);
    MPI_Send(sendbuf,nx_local*(nz_local+2),MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++)
	lattice[i][ny_local][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (k = 0; k <= nz_local+1; k++) {
	lattice[i][1][k] = lattice[i][ny_local+1][k];
	lattice[i][ny_local][k] = lattice[i][0][k];
      }
  }

  // send down, then up
  // if only 1 proc in z, copy data

  if (procdown != me) {
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	sendbuf[m++] = lattice[i][j][nz_local+1];

    MPI_Irecv(recvbuf,nx_local*ny_local,MPI_INT,procup,0,world,&request);
    MPI_Send(sendbuf,nx_local*ny_local,MPI_INT,procdown,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	lattice[i][j][1] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	sendbuf[m++] = lattice[i][j][0];

    MPI_Irecv(recvbuf,nx_local*ny_local,MPI_INT,procdown,0,world,&request);
    MPI_Send(sendbuf,nx_local*ny_local,MPI_INT,procup,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	lattice[i][j][nz_local] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++) {
	lattice[i][j][1] = lattice[i][j][nz_local+1];
	lattice[i][j][nz_local] = lattice[i][j][0];
      }
  }

}

/* ----------------------------------------------------------------------
   send back ghost values for entire sub-domain owned by this proc
   dellocal > 1
------------------------------------------------------------------------- */

void CommLattice3d::reverse_all_multilayer(int ***lattice)
{
  int i,j,k,m,nsend;
  MPI_Request request;
  MPI_Status status;


  // send west, then east
  // if only 1 proc in x, copy data

  if (procwest != me) {
    nsend = (ny_local+2*dellocal)*(nz_local+2*dellocal)*dellocal;
    m = 0;
    for (i = 1; i <= dellocal; i++)
      for (j = 1-dellocal; j <= ny_local+dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  sendbuf[m++] = lattice[nx_local+i][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,proceast,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procwest,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= dellocal; i++)
      for (j = 1-dellocal; j <= ny_local+dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  lattice[i][j][k] = recvbuf[m++];


    m = 0;
    for (i = 1; i <= dellocal; i++)
      for (j = 1-dellocal; j <= ny_local+dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  sendbuf[m++] = lattice[i-dellocal][j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procwest,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,proceast,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= dellocal; i++)
      for (j = 1-dellocal; j <= ny_local+dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  lattice[nx_local-dellocal+i][j][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= dellocal; i++)
      for (j = 1-dellocal; j <= ny_local+dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++) {
	  lattice[i][j][k] = lattice[nx_local+i][j][k];
	  lattice[nx_local-dellocal+i][j][k] = lattice[i-dellocal][j][k];
	}
  }

  // send south, then north
  // if only 1 proc in y, copy data

  if (procsouth != me) {
    nsend = nx_local*(nz_local+2*dellocal)*dellocal;
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  sendbuf[m++] = lattice[i][ny_local+j][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procnorth,0,world,
	      &request);
    MPI_Send(sendbuf,nsend,MPI_INT,procsouth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  lattice[i][j][k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  sendbuf[m++] = lattice[i][j-dellocal][k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procsouth,0,world,
	      &request);
    MPI_Send(sendbuf,nsend,MPI_INT,procnorth,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++)
	  lattice[i][ny_local-dellocal+j][k] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= dellocal; j++)
	for (k = 1-dellocal; k <= nz_local+dellocal; k++) {
	  lattice[i][j][k] = lattice[i][ny_local+j][k];
	  lattice[i][ny_local-dellocal+j][k] = lattice[i][j-dellocal][k];
	}
  }

  // send down, then up
  // if only 1 proc in z, copy data

  if (procdown != me) {
    nsend = nx_local*ny_local*dellocal;
    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= dellocal; k++)
	sendbuf[m++] = lattice[i][j][nz_local+k];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procup,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procdown,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= dellocal; k++)
	  lattice[i][j][k] = recvbuf[m++];

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= dellocal; k++)
	  sendbuf[m++] = lattice[i][j][k-dellocal];

    MPI_Irecv(recvbuf,nsend,MPI_INT,procdown,0,world,&request);
    MPI_Send(sendbuf,nsend,MPI_INT,procup,0,world);
    MPI_Wait(&request,&status);

    m = 0;
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= dellocal; k++)
	  lattice[i][j][nz_local-dellocal+k] = recvbuf[m++];

  } else {
    for (i = 1; i <= nx_local; i++)
      for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= dellocal; k++) {
	lattice[i][j][k] = lattice[i][j][nz_local+k];
	lattice[i][j][nz_local-dellocal+k] = lattice[i][j][k-dellocal];
	}
  }
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice3d::allocate_swap(const int idim, const int jdim)
{
  memory->create_2d_T_array(swapinfo,idim,jdim,"commlattice:swapinfo");
  memory->create_2d_T_array(reverseinfo,idim,jdim,"commlattice:reverseinfo");
  maxsector = idim;
  maxswap = jdim;
}
