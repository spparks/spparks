/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_energy.h"
#include "app_lattice.h"
#include "comm_lattice.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEnergy::DiagEnergy(SPK *spk, int narg, char **arg) : Diag(spk,narg,arg)
{

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style energy output file");
	}
      } else {
	error->all("Illegal diag_style energy command");
      } 
    } else {
      //      error->all("Illegal diag_style energy command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagEnergy::~DiagEnergy()
{
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::init(double time)
{
  if (app->appclass != App::LATTICE) error->all("diag_style incompatible with app_style");

  applattice = (AppLattice *) app;

  nlocal = applattice->nlocal;

  applattice->comm->all();
  
  double etmp = 0.0;
  for (int i = 0; i < nlocal; i++)
    etmp += applattice->site_energy(i);
  
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagEnergy::compute(double time, int done)
{
  double etmp;

  if (check_time(time, done))  {

    applattice->comm->all();

    etmp = 0.0;
    for (int i = 0; i < nlocal; i++) 
      etmp += applattice->site_energy(i);

    
    MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  }
  
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats(char *strtmp) {
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats_header(char *strtmp) {
  sprintf(strtmp," %10s","Energy");
}
