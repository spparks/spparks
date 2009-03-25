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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "diag.h"
#include "app.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Diag::Diag(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  // stats_flag = 1, provide stats, compute interval controlled by Output
  // stats_flag = 0, do not provide stats to Output

  stats_flag = 1;
  diag_delta = 0.0;
  diag_logfreq = 0;
  diag_delay = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if  (strcmp(arg[iarg],"stats") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"yes") == 0) stats_flag = 1;
	else if (strcmp(arg[iarg],"no") == 0) stats_flag = 0;
	else error->all("Illegal diag_style command");
      } else error->all("Illegal diag_style command");
    } else if (strcmp(arg[iarg],"delt") == 0) {
      iarg++;
      if (iarg < narg) {
	diag_delta = atof(arg[iarg]);
	if (diag_delta <= 0.0) error->all("Illegal diag_style command");
      } else error->all("Illegal diag_style command");
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      diag_logfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	diag_nrepeat = atoi(arg[iarg]);
	iarg++;
	diag_scale = atof(arg[iarg]);
	if (diag_scale <= 0.0) error->all("Illegal diag_style command");
	if (diag_nrepeat*diag_delta > diag_scale)
	  error->all("Illegal diag_style command");
      } else error->all("Illegal diag_style command");
    } else if (strcmp(arg[iarg],"delay") == 0) {
      iarg++;
      if (iarg < narg) {
	diag_delay = atof(arg[iarg]);
      } else error->all("Illegal diag_style command");
    } else break;
    iarg++;
  }

  iarg_child = iarg;

  if (stats_flag && diag_logfreq) error->all("Illegal diag_style command");
  if (stats_flag && diag_delta > 0.0) error->all("Illegal diag_style command");
  if (stats_flag && diag_delay > 0.0) error->all("Illegal diag_style command");
}

/* ---------------------------------------------------------------------- */

Diag::~Diag()
{
  delete [] style;
}
