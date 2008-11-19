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

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  // Default set stats_flag, so that stats are provided to Output 
  // and output interval is controlled by Output
  //
  // stats_flag = 0, do not provide stats to Output
  // stats_flag = 1, provide stats, compute interval controlled by Output
  //
  // If stats_flag = 1, then require diag_delta = 0.0; 

  stats_flag = 1;
  diag_delta = 0.0;
  diag_ilogfreq = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"delt") == 0) {
      iarg++;
      if (iarg < narg) {
	diag_delta = atof(arg[iarg]);
      } else error->all("Illegal diag_style command");
    } else if  (strcmp(arg[iarg],"stats") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"yes") == 0) stats_flag = 1;
	else if (strcmp(arg[iarg],"no") == 0) stats_flag = 0;
	else error->all("Illegal diag_style command");
      } else error->all("Illegal diag_style command");
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      diag_ilogfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	diag_nrepeat = atoi(arg[iarg]);
	iarg++;
	diag_scale = atof(arg[iarg]);
      } else error->all("Illegal diag_style command");
    } else break;
    iarg++;
  }

  iarg_child = iarg;

  if (diag_delta < 0.0) error->all("Illegal diag_style command");
  if (diag_ilogfreq && diag_delta <= 0.0) 
    error->all("Illegal diag_style command");
  if (stats_flag && diag_delta > 0.0)
    error->all("Illegal diag_style command");
}

/* ---------------------------------------------------------------------- */

Diag::~Diag()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

int Diag::check_time(double time, int done)
{
  int iflag = 0;
  if ((diag_delta > 0.0 && time >= diag_time) || done) iflag = 1;

  if ((diag_delta > 0.0 && time >= diag_time)) {
    if (diag_ilogfreq == 0) {
      // Ensure that new diag_time exceeds time
      while (time >= diag_time) {
	diag_time += diag_delta;
      }
      
    } else if (diag_ilogfreq == 1) {
      // Ensure that new diag_time exceeds time
      while (time >= diag_time) {
	diag_time += diag_delta;
	diag_irepeat++;
	if (diag_irepeat == diag_nrepeat) {
	  diag_delta *= diag_scale;
	  diag_time = diag_t0+diag_delta;
	  diag_irepeat = 0;
	}
      }
    }
    
  }
  return iflag;
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

