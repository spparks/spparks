/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_grain_3d.h"
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

CommGrain3D::CommGrain3D(class SPK *spk) : 
  CommGrain(spk),
  sendbuf(NULL),
  recvbuf(NULL),
  swapinfo(NULL)
{}

/* ---------------------------------------------------------------------- */

CommGrain3D::~CommGrain3D() 
{
  free_swap();
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommGrain3D::setup(const int nx, const int ny, const int nz, 
			 const int nxh, const int nyh, const int nzh, 
			 const int procwest, const int proceast, 
			 const int procsouth, const int procnorth,
			 const int procdown, const int procup)
{
  // initialize swap parameters and allocate memory
  // This 3D scheme for assigning swap parameters is more automated than
  // the one used for 2D. This makes it easy to modify for different applications
  // It also makes it possible to extend to higher-dimensional spaces!! 
  SwapInfo* swap;
  int iswap,isector;
  int i,j,k,ii,jj,kk,irecv,isend;
  int procarray[3][2];
  SwapInfo swapguide[2][2];

  // isector = 0,1,2,3 refer to lower SW,NW,SE,NE quadrants, respectively
  // isector = 4,5,6,7 refer to uper SW,NW,SE,NE quadrants, respectively
  // iswap = 0 refers to the up/down swap
  // iswap = 1 refers to the north/south swap
  // iswap = 2 refers to the east/west swap

  nsector = 8;
  nswap = 3;
  allocate_swap(nsector,nswap);

  // Assign processors to different swap directions as [iswap=0,1,2][dir=0,1]
  procarray[0][0] = procdown;
  procarray[0][1] = procup;
  procarray[1][0] = procsouth;
  procarray[1][1] = procnorth;
  procarray[2][0] = procwest;
  procarray[2][1] = proceast;

  // swapguide[0][0] stores values for lower side, not communicating
  swapguide[0][0].recvix = 0;
  swapguide[0][0].recviy = 0;
  swapguide[0][0].recviz = 0;
  swapguide[0][0].sendix = 0;
  swapguide[0][0].sendiy = 0;
  swapguide[0][0].sendiz = 0;
  swapguide[0][0].numx = nxh+1;
  swapguide[0][0].numy = nyh+1;
  swapguide[0][0].numz = nzh+1;
  swapguide[0][0].recvproc = -1;
  swapguide[0][0].sendproc = -1;
  swapguide[0][0].copyixa = nx;
  swapguide[0][0].copyixb = 0;
  swapguide[0][0].copyiya = ny;
  swapguide[0][0].copyiyb = 0;

  // swapguide[0][1] stores values for lower side, communicating
  swapguide[0][1].recvix = 0;
  swapguide[0][1].recviy = 0;
  swapguide[0][1].recviz = 0;
  swapguide[0][1].sendix = nx;
  swapguide[0][1].sendiy = ny;
  swapguide[0][1].sendiz = nz;
  swapguide[0][1].numx = 1;
  swapguide[0][1].numy = 1;
  swapguide[0][1].numz = 1;
  swapguide[0][1].recvproc = 0;
  swapguide[0][1].sendproc = 1;
  swapguide[0][1].copyixa = nx;
  swapguide[0][1].copyixb = 0;
  swapguide[0][1].copyiya = ny;
  swapguide[0][1].copyiyb = 0;

  // swapguide[1][0] stores values for upper side, not communicating
  swapguide[1][0].recvix = nxh-1;
  swapguide[1][0].recviy = nyh-1;
  swapguide[1][0].recviz = nzh-1;
  swapguide[1][0].sendix = nxh-1;
  swapguide[1][0].sendiy = nyh-1;
  swapguide[1][0].sendiz = nzh-1;
  swapguide[1][0].numx = nx-nxh+3;
  swapguide[1][0].numy = ny-nyh+3;
  swapguide[1][0].numz = nz-nzh+3;
  swapguide[1][0].recvproc = -1;
  swapguide[1][0].sendproc = -1;
  swapguide[1][0].copyixa = 1;
  swapguide[1][0].copyixb = nx+1;
  swapguide[1][0].copyiya = 1;
  swapguide[1][0].copyiyb = ny+1;

  // swapguide[1][1] stores values for upper side, communicating
  swapguide[1][1].recvix = nx+1;
  swapguide[1][1].recviy = ny+1;
  swapguide[1][1].recviz = nz+1;
  swapguide[1][1].sendix = 1;
  swapguide[1][1].sendiy = 1;
  swapguide[1][1].sendiz = 1;
  swapguide[1][1].numx = 1;
  swapguide[1][1].numy = 1;
  swapguide[1][1].numz = 1;
  swapguide[1][1].recvproc = 1;
  swapguide[1][1].sendproc = 0;
  swapguide[1][1].copyixa = 1;
  swapguide[1][1].copyixb = nx+1;
  swapguide[1][1].copyiya = 1;
  swapguide[1][1].copyiyb = ny+1;

  // Loop over both sides in each direction, and all swap directions
  // Use appropriate guide values in each case
  isector = 0;
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      for (int k=0;k<2;k++) {
	for (iswap=0;iswap<nswap;iswap++) {
	  swap = &swapinfo[isector][iswap];

	  // Identify which direction is communicating	  
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
  }
  
  maxrecv = 0;
  maxsend = 0;
  for (isector=0;isector<nsector;isector++) {
    for (iswap=0;iswap<nswap;iswap++) {
      swap = &swapinfo[isector][iswap];
      maxrecv=MAX(maxrecv,swap->numx*swap->numy*swap->numz);
      maxsend=MAX(maxsend,swap->numx*swap->numy*swap->numz);
    }
  }

  recvbuf = (int*) memory->smalloc(maxrecv*sizeof(int),"commgrain:recvbuf");
  sendbuf = (int*) memory->smalloc(maxsend*sizeof(int),"commgrain:sendbuf");


}

/* ----------------------------------------------------------------------
   communication of boundary values for one sector
------------------------------------------------------------------------- */

void CommGrain3D::communicate(int*** lattice, const int isector) {
  int iswap,ix,iy,iz,ixa,ixb,iya,iyb,ii,i,j,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;

  // Each swap direction is a little different
  // Need to code each swap explicitly

  // First swap is in the z direction
  iswap = 0;
  swap = &swapinfo[isector][iswap];
  recvnum = swap->numx*swap->numy;
  sendnum = recvnum;
  ixa = swap->copyixa;
  ixb = swap->copyixb;
  iya = swap->copyiya;
  iyb = swap->copyiyb;
  iz = swap->sendiz;

  // Copy y-edge to ghost edge

  iy = swap->sendiy;
  for (i=0;i<swap->numy;i++) {
    lattice[ixb][iy][iz] = lattice[ixa][iy][iz];
    iy++;
  }

  // Copy x-edge to ghost edge

  ix = swap->sendix;
  for (i=0;i<swap->numx;i++) {
    lattice[ix][iyb][iz] = lattice[ix][iya][iz];
    ix++;
  }

  // Copy corner to ghost corner; must do this last to avoid overwrite

  lattice[ixb][iyb][iz] = lattice[ixa][iya][iz];

  // Use MPI to pass data only with other processes

  if (swap->sendproc != universe->me) {

    // Pack the send buffer

    ii = 0;
    ix = swap->sendix;
    for (i=0;i<swap->numx;i++) {
      iy = swap->sendiy;
      for (j=0;j<swap->numy;j++) {
	sendbuf[ii] = lattice[ix][iy][iz];
	iy++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    iz = swap->recviz;
    ix = swap->recvix;
    for (i=0;i<swap->numx;i++) {
      iy = swap->recviy;
      for (j=0;j<swap->numy;j++) {
	lattice[ix][iy][iz] = recvbuf[ii]; 
	iy++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    iz = swap->recviz;
    ix = swap->recvix;
    for (i=0;i<swap->numx;i++) {
      iy = swap->recviy;
      for (j=0;j<swap->numy;j++) {
	lattice[ix][iy][iz] = lattice[ix][iy][swap->sendiz]; 
	iy++;
      }
      ix++;
    }
  }

  iz = swap->recviz;

  // Copy ghost y-edge to edge

  iy = swap->recviy;
  for (i=0;i<swap->numy;i++) {
    lattice[ixa][iy][iz] = lattice[ixb][iy][iz];
    iy++;
  }

  // Copy ghost x-edge to edge

  ix = swap->recvix;
  for (i=0;i<swap->numx;i++) {
    lattice[ix][iya][iz] = lattice[ix][iyb][iz];
    ix++;
  }

  // Copy ghost corner to corner; must do this last to avoid overwrite

  lattice[ixa][iya][iz] = lattice[ixb][iyb][iz];

  // Second swap is in the y direction
  iswap = 1;
  swap = &swapinfo[isector][iswap];
  recvnum = swap->numx*swap->numz;
  sendnum = recvnum;
  iy = swap->sendiy;
  
  // Copy edge to ghost edge

  ixa = swap->copyixa;
  ixb = swap->copyixb;
  iz = swap->sendiz;
  for (i=0;i<swap->numz;i++) {
    lattice[ixb][iy][iz] = lattice[ixa][iy][iz];
    iz++;
  }

  // Use MPI to pass data only with other processes

  if (swap->sendproc != universe->me) {

    // Pack the send buffer

    ii = 0;
    ix = swap->sendix;
    for (i=0;i<swap->numx;i++) {
      iz = swap->sendiz;
      for (j=0;j<swap->numz;j++) {
	sendbuf[ii] = lattice[ix][iy][iz];
	iz++;
	ii++;
      }
      ix++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    iy = swap->recviy;
    ix = swap->recvix;
    for (i=0;i<swap->numx;i++) {
      iz = swap->recviz;
      for (j=0;j<swap->numz;j++) {
	lattice[ix][iy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      ix++;
    }

  // Use direct access to pass data to self

  } else { 

    iy = swap->recviy;
    ix = swap->recvix;
    for (i=0;i<swap->numx;i++) {
      iz = swap->recviz;
      for (j=0;j<swap->numz;j++) {
	lattice[ix][iy][iz] = lattice[ix][swap->sendiy][iz]; 
	iz++;
      }
      ix++;
    }
  }

  iy = swap->recviy;

  // Copy ghost edge to edge

  ixa = swap->copyixa;
  ixb = swap->copyixb;
  iz = swap->recviz;
  for (i=0;i<swap->numz;i++) {
    lattice[ixa][iy][iz] = lattice[ixb][iy][iz];
    iz++;
  }

  // Third swap is in the x direction

  iswap = 2;
  swap = &swapinfo[isector][iswap];
  recvnum = swap->numy*swap->numz;
  sendnum = recvnum;
  ix = swap->sendix;

  // Use MPI to pass data only with other processes

  if (swap->sendproc != universe->me) {

    // Pack the send buffer

    ii = 0;
    iy = swap->sendiy;
    for (i=0;i<swap->numy;i++) {
      iz = swap->sendiz;
      for (j=0;j<swap->numz;j++) {
	sendbuf[ii] = lattice[ix][iy][iz];
	iz++;
	ii++;
      }
      iy++;
    }

    MPI_Irecv(recvbuf,recvnum,MPI_INT,
	      swap->recvproc,0,world,&request);
    MPI_Send(sendbuf,sendnum,MPI_INT,
	     swap->sendproc,0,world);
    MPI_Wait(&request,&status);
    
    // Unpack the receive buffer

    ii = 0;
    ix = swap->recvix;
    iy = swap->recviy;
    for (i=0;i<swap->numy;i++) {
      iz = swap->recviz;
      for (j=0;j<swap->numz;j++) {
	lattice[ix][iy][iz] = recvbuf[ii]; 
	iz++;
	ii++;
      }
      iy++;
    }

  // Use direct access to pass data to self

  } else { 

    ix = swap->recvix;
    iy = swap->recviy;
    for (i=0;i<swap->numy;i++) {
      iz = swap->recviz;
      for (j=0;j<swap->numz;j++) {
	lattice[ix][iy][iz] = lattice[swap->sendix][iy][iz]; 
	iz++;
      }
      iy++;
    }
  }

}

/* ----------------------------------------------------------------------
   deallocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommGrain3D::free_swap()
{
  memory->destroy_2d_T_array(swapinfo);
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommGrain3D::allocate_swap(const int idim, const int jdim)
{
  // Allocate a 2D array of SwapInfo objects

  memory->create_2d_T_array(swapinfo,idim,jdim,"commgrain3d:swapinfo");
  maxsector = idim;
  maxswap = jdim;
}
