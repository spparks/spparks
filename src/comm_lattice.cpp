/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "comm_lattice.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

CommLattice::CommLattice(class SPK *spk) : SysPtr(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  //swapinfo = NULL;
  //reverseinfo = NULL;
}

/* ---------------------------------------------------------------------- */

CommLattice::~CommLattice()
{
  free_swap();
  //memory->sfree(sendbuf);
  //memory->sfree(recvbuf);
}

/* ---------------------------------------------------------------------- */

void CommLattice::init(const int nlocal_in, 
		       const int procwest_in, const int proceast_in, 
		       const int procsouth_in, const int procnorth_in,
		       const int delghost_in, const int dellocal_in)
{
  dimension = 2;
  nlocal = nlocal_in;

  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;

  delghost = delghost_in;
  dellocal = dellocal_in;
}

/* ---------------------------------------------------------------------- */

void CommLattice::init(const int nlocal_in, 
		       const int procwest_in, const int proceast_in, 
		       const int procsouth_in, const int procnorth_in,
		       const int procdown_in, const int procup_in,
		       const int delghost_in, const int dellocal_in)
{
  dimension = 3;
  nlocal = nlocal_in;

  procwest = procwest_in;
  proceast = proceast_in;
  procsouth = procsouth_in;
  procnorth = procnorth_in;
  procdown = procdown_in;
  procup = procup_in;

  delghost = delghost_in;
  dellocal = dellocal_in;
}

/* ----------------------------------------------------------------------
   setup up forward communication for each sector
------------------------------------------------------------------------- */

void CommLattice::setup_swapinfo()
{
}

/* ----------------------------------------------------------------------
   setup up reverse communication for each sector
------------------------------------------------------------------------- */

void CommLattice::setup_reverseinfo()
{
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   delghost > 1
------------------------------------------------------------------------- */

void CommLattice::sector_multilayer(int* lattice, const int isector)
{
}

/* ----------------------------------------------------------------------
   communicate ghost values for one sector (quadrant)
   delghost > 1
   this version erases all local sites that are read
------------------------------------------------------------------------- */

void CommLattice::sector_multilayer_destroy(int* lattice, const int isector)
{
}

/* ----------------------------------------------------------------------
   send back ghost values for one sector (quadrant)
   it chooses the only variant: reverse_sector_multilayer
------------------------------------------------------------------------- */

void CommLattice::reverse_sector(int* lattice, const int isector) 
{
}

/* ----------------------------------------------------------------------
   send back ghost values for one sector (quadrant)
   delghost > 1
------------------------------------------------------------------------- */

void CommLattice::reverse_sector_multilayer(int* lattice, const int isector)
{
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
   it chooses the only variant: all_multilayer
------------------------------------------------------------------------- */

void CommLattice::all(int* lattice)
{
}

/* ----------------------------------------------------------------------
   update ghost values for entire sub-domain owned by this proc
   delghost > 1
------------------------------------------------------------------------- */

void CommLattice::all_multilayer(int *lattice)
{
}

/* ----------------------------------------------------------------------
   deallocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice::free_swap()
{
}

/* ----------------------------------------------------------------------
   allocate arrays for swap parameters
------------------------------------------------------------------------- */

void CommLattice::allocate_swap(const int idim, const int jdim)
{
}
