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
}

/* ---------------------------------------------------------------------- */

CommLattice3d::~CommLattice3d()
{
  free_swap();
  memory->sfree(sendbuf);
  memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommLattice3d::init(const int nx_local_in, const int ny_local_in,
			 const int nz_local_in,
			 const int procwest_in, const int proceast_in, 
			 const int procsouth_in, const int procnorth_in,
			 const int procdown_in, const int procup_in)
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

  int nx_half = nx_local/2 + 1;
  int ny_half = ny_local/2 + 1;
  int nz_half = nz_local/2 + 1;

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
  swapguide[0][0].numx = nx_half+1;
  swapguide[0][0].numy = ny_half+1;
  swapguide[0][0].numz = nz_half+1;
  swapguide[0][0].recvproc = -1;
  swapguide[0][0].sendproc = -1;
  swapguide[0][0].copyixa = nx_local;
  swapguide[0][0].copyixb = 0;
  swapguide[0][0].copyiya = ny_local;
  swapguide[0][0].copyiyb = 0;

  // swapguide[0][1] stores values for lower side, communicating

  swapguide[0][1].recvix = 0;
  swapguide[0][1].recviy = 0;
  swapguide[0][1].recviz = 0;
  swapguide[0][1].sendix = nx_local;
  swapguide[0][1].sendiy = ny_local;
  swapguide[0][1].sendiz = nz_local;
  swapguide[0][1].numx = 1;
  swapguide[0][1].numy = 1;
  swapguide[0][1].numz = 1;
  swapguide[0][1].recvproc = 0;
  swapguide[0][1].sendproc = 1;
  swapguide[0][1].copyixa = nx_local;
  swapguide[0][1].copyixb = 0;
  swapguide[0][1].copyiya = ny_local;
  swapguide[0][1].copyiyb = 0;

  // swapguide[1][0] stores values for upper side, not communicating

  swapguide[1][0].recvix = nx_half-1;
  swapguide[1][0].recviy = ny_half-1;
  swapguide[1][0].recviz = nz_half-1;
  swapguide[1][0].sendix = nx_half-1;
  swapguide[1][0].sendiy = ny_half-1;
  swapguide[1][0].sendiz = nz_half-1;
  swapguide[1][0].numx = nx_local-nx_half+3;
  swapguide[1][0].numy = ny_local-ny_half+3;
  swapguide[1][0].numz = nz_local-nz_half+3;
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
  swapguide[1][1].numx = 1;
  swapguide[1][1].numy = 1;
  swapguide[1][1].numz = 1;
  swapguide[1][1].recvproc = 1;
  swapguide[1][1].sendproc = 0;
  swapguide[1][1].copyixa = 1;
  swapguide[1][1].copyixb = nx_local+1;
  swapguide[1][1].copyiya = 1;
  swapguide[1][1].copyiyb = ny_local+1;

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

  maxbuf = (nx_local+2) * (ny_local*2);
  maxbuf = MAX(maxbuf,(nx_local+2) * (nz_local*2));
  maxbuf = MAX(maxbuf,(ny_local+2) * (nz_local*2));
  recvbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:recvbuf");
  sendbuf = (int*) memory->smalloc(maxbuf*sizeof(int),"commlattice:sendbuf");
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
------------------------------------------------------------------------- */

void CommLattice3d::sector(int*** lattice, const int isector) {
  int iswap,ii,ix,iy,iz,i,j,recvnum,sendnum;
  MPI_Request request;
  MPI_Status status;
  SwapInfo* swap;

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
   update ghost values for entire sub-domain owned by this proc
------------------------------------------------------------------------- */

void CommLattice3d::all(int ***lattice)
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
   deallocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice3d::free_swap()
{
  memory->destroy_2d_T_array(swapinfo);
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice3d::allocate_swap(const int idim, const int jdim)
{
  memory->create_2d_T_array(swapinfo,idim,jdim,"commlattice:swapinfo");
  maxsector = idim;
  maxswap = jdim;
}
