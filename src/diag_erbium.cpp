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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_erbium.h"
#include "app.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagErbium::DiagErbium(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
}

/* ---------------------------------------------------------------------- */

void DiagErbium::init(double time)
{
  if (app->appclass != App::LATTICE)
    error->all("Diag style incompatible with app style");

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

void DiagErbium::compute(double time, int iflag, int done)
{
  double etmp;

  if (stats_flag == 0) {
    iflag = check_time(time, done);
  }

  if (iflag || done) {
    applattice->comm->all();

    etmp = 0.0;
    for (int i = 0; i < nlocal; i++) 
      etmp += applattice->site_energy(i);
    
    MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Energy");
}
