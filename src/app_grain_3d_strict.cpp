/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_grain_3d.h"
#include "app_grain_3d_strict.h"
#include "comm_grain_3d.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace SPPARKS;
using namespace std;

/* ---------------------------------------------------------------------- */

AppGrain3DStrict::AppGrain3DStrict(SPK *spk, int narg, char **arg) : AppGrain3D(spk,narg,arg)
{
  int i,j,k,ii,jj,kk;
  int isite;

  // Unfortunately, this error will never be reached, 
  // because it will be preceeded by error in AppGrain3D constructor
  if (narg != 6) error->all("Invalid app_style grain_3d_strict command");

  // init ghost spins (take this out later)

  for (i = 0; i <= nx_local+1; i++) 
    for (j = 0; j <= ny_local+1; j++) 
      for (k = 0; k <= nz_local+1; k++) 
      lattice[i][j][k] = 0;
  
  // initialize spins
  
  random = new RandomPark(seed);
  for (i = 0; i < 100; i++) random->uniform();
  
  // loop over global list
  // so that assigment is indpendent of parallel decompostion
  for (i = 1; i <= nx_global; i++) {
    ii = i - nx_offset;
    if (ii >= 1 && ii <= nx_local) { 
      for (j = 1; j <= ny_global; j++) {
	jj = j - ny_offset;
	if (jj >= 1 && jj <= ny_local) { 
	  for (k = 1; k <= nz_global; k++) {
	    kk = k - nz_offset;
	    if (kk >= 1 && kk <= nz_local) { 
	      lattice[ii][jj][kk] = random->irandom(nspins);
	    } else {
	      random->irandom(nspins);
	    }
	  }
	} else {
	  for (k = 1; k <= nz_global; k++) {
	    random->irandom(nspins);
	  }
	}
      }
    } else {
      for (j = 1; j <= ny_global; j++) {
	for (k = 1; k <= nz_global; k++) {
	  random->irandom(nspins);
	}
      }
    }
  }

  // Set up array of random number generators
  // Each randomPark object uses only 4 bytes 
  memory->create_3d_T_array(ranlat,nx_local+2,ny_local+2,nz_local+2,
			      "grain3d:ranlat");

  // construct random number generators

  for (i = 1; i <= nx_local; i++) {
    for (j = 1; j <= ny_local; j++) {
      for (k = 1; k <= nz_local; k++) {
	isite = ((i+nx_offset)*ny_global + j + ny_offset)*nz_global+ k + nz_offset;
	ranlat[i][j][k].init(seed+isite);
      }
    }
  }
  
}

/* ---------------------------------------------------------------------- */

AppGrain3DStrict::~AppGrain3DStrict()
{
  memory->destroy_3d_T_array(ranlat);
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppGrain3DStrict::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  nsweep = atoi(arg[0]);
  
  // error check
  
  //if (solve == NULL) error->all("No solver class defined");
  
  // init classes used by this app
  
  int i;
  double *tmp;
  
  init();
  timer->init();
  
  // perform the run
  
  iterate();
  
  // final statistics
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on sweep solver; event order strictly matches serial case
 ------------------------------------------------------------------------- */

void AppGrain3DStrict::iterate()
{
  timer->barrier_start(TIME_LOOP);
  
  for (int i = 0; i < nsweep; i++) {
    ntimestep++;

    // Loop over the lattice colors a la checkerboarding for 2D Ising model
    for (int icolor = 0; icolor < nsector; icolor++) {
      // Loop over the sectors
      for (int isector = 0; isector < nsector; isector++) {
	timer->stamp();
	comm->communicate(lattice,isector);
	timer->stamp(TIME_COMM);

	timer->stamp();
	sweep(icolor,isector);
	timer->stamp(TIME_SOLVE);
      }
    }

    if (ntimestep == stats_next) {
      timer->stamp();
      stats();
      stats_next += nstats;
      timer->stamp(TIME_OUTPUT);
    }
    if (ntimestep == dump_next) {
      timer->stamp();
      dump();
      dump_next += ndump;
      timer->stamp(TIME_OUTPUT);
    }
  }
  
  timer->barrier_stop(TIME_LOOP);
}

/* -----------------------------------------------------------------------------
   perform spin flips for one color in one sector; event order strictly matches serial case
 ------------------------------------------------------------------------------- */
void AppGrain3DStrict::sweep(int icolor, int isector)
{
  int i,j,k,iold,inew,nold,nnew,i0,j0,k0;

  int xlo = quad[isector].xlo;
  int xhi = quad[isector].xhi;
  int ylo = quad[isector].ylo;
  int yhi = quad[isector].yhi;
  int zlo = quad[isector].zlo;
  int zhi = quad[isector].zhi;

  i0 = (icolor/4 + nx_offset + xlo) % 2;
  j0 = (icolor/2 + ny_offset + ylo) % 2;
  k0 = (icolor/1 + nz_offset + zlo) % 2;

  for (i = xlo+i0; i <= xhi; i+=2) {
    for (j = ylo+j0; j <= yhi; j+=2) {
      for (k = zlo+k0; k <= zhi; k+=2) {

	iold = lattice[i][j][k];
	nold = 0;

	if (iold == lattice[i-1][j-1][k-1]) nold++;
	if (iold == lattice[i-1][j-1][k]) nold++;
	if (iold == lattice[i-1][j-1][k+1]) nold++;

	if (iold == lattice[i-1][j][k-1]) nold++;
	if (iold == lattice[i-1][j][k]) nold++;
	if (iold == lattice[i-1][j][k+1]) nold++;

	if (iold == lattice[i-1][j+1][k-1]) nold++;
	if (iold == lattice[i-1][j+1][k]) nold++;
	if (iold == lattice[i-1][j+1][k+1]) nold++;

	if (iold == lattice[i][j-1][k-1]) nold++;
	if (iold == lattice[i][j-1][k]) nold++;
	if (iold == lattice[i][j-1][k+1]) nold++;

	if (iold == lattice[i][j][k-1]) nold++;
	if (iold == lattice[i][j][k+1]) nold++;

	if (iold == lattice[i][j+1][k-1]) nold++;
	if (iold == lattice[i][j+1][k]) nold++;
	if (iold == lattice[i][j+1][k+1]) nold++;

	if (iold == lattice[i+1][j-1][k-1]) nold++;
	if (iold == lattice[i+1][j-1][k]) nold++;
	if (iold == lattice[i+1][j-1][k+1]) nold++;

	if (iold == lattice[i+1][j][k-1]) nold++;
	if (iold == lattice[i+1][j][k]) nold++;
	if (iold == lattice[i+1][j][k+1]) nold++;

	if (iold == lattice[i+1][j+1][k-1]) nold++;
	if (iold == lattice[i+1][j+1][k]) nold++;
	if (iold == lattice[i+1][j+1][k+1]) nold++;

	inew = ranlat[i][j][k].irandom(nspins);
	nnew = 0;

	if (inew == lattice[i-1][j-1][k-1]) nnew++;
	if (inew == lattice[i-1][j-1][k]) nnew++;
	if (inew == lattice[i-1][j-1][k+1]) nnew++;

	if (inew == lattice[i-1][j][k-1]) nnew++;
	if (inew == lattice[i-1][j][k]) nnew++;
	if (inew == lattice[i-1][j][k+1]) nnew++;

	if (inew == lattice[i-1][j+1][k-1]) nnew++;
	if (inew == lattice[i-1][j+1][k]) nnew++;
	if (inew == lattice[i-1][j+1][k+1]) nnew++;

	if (inew == lattice[i][j-1][k-1]) nnew++;
	if (inew == lattice[i][j-1][k]) nnew++;
	if (inew == lattice[i][j-1][k+1]) nnew++;

	if (inew == lattice[i][j][k-1]) nnew++;
	if (inew == lattice[i][j][k+1]) nnew++;

	if (inew == lattice[i][j+1][k-1]) nnew++;
	if (inew == lattice[i][j+1][k]) nnew++;
	if (inew == lattice[i][j+1][k+1]) nnew++;

	if (inew == lattice[i+1][j-1][k-1]) nnew++;
	if (inew == lattice[i+1][j-1][k]) nnew++;
	if (inew == lattice[i+1][j-1][k+1]) nnew++;

	if (inew == lattice[i+1][j][k-1]) nnew++;
	if (inew == lattice[i+1][j][k]) nnew++;
	if (inew == lattice[i+1][j][k+1]) nnew++;

	if (inew == lattice[i+1][j+1][k-1]) nnew++;
	if (inew == lattice[i+1][j+1][k]) nnew++;
	if (inew == lattice[i+1][j+1][k+1]) nnew++;

	if (nold <= nnew) lattice[i][j][k] = inew;
      }
    }
  }
}

