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
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Diag::Diag(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  diag_ilogfreq = 0;

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  if (narg < 2) error->all("Illegal diag_style command");
  diag_delta = atof(arg[1]);
  if (diag_delta <= 0.0) error->all("Illegal diag_style command");

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      diag_ilogfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	diag_nrepeat = atoi(arg[iarg]);
	iarg++;
	diag_scale = atof(arg[iarg]);
      } else {
	error->all("Illegal diag_style command");
      }
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Diag::~Diag()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

int Diag::check_time(double time, int done)
{
  if ((diag_delta > 0.0 && time >= diag_time) || done) {
    if (diag_ilogfreq == 0) {
      diag_time += diag_delta;
    } else if (diag_ilogfreq == 1) {
      diag_time += diag_delta;
      diag_irepeat++;
      if (diag_irepeat == diag_nrepeat) {
	diag_delta *= diag_scale;
	diag_time = diag_t0+diag_delta;
	diag_irepeat = 0;
      }
    }
    return 1;
  } else {
    return 0;
  }
}

/* ---------------------------------------------------------------------- */

void Diag::setup_time(double time)
{
  if (diag_ilogfreq == 0) {
    diag_time = time + diag_delta;
  } else if (diag_ilogfreq == 1) {
    diag_time = time + diag_delta;
    diag_t0 = time;
    diag_irepeat = 0;
  }
}

