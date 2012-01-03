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
  delta = 0.0;
  logfreq = 0;
  delay = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stats") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"yes") == 0) stats_flag = 1;
	else if (strcmp(arg[iarg],"no") == 0) stats_flag = 0;
	else error->all(FLERR,"Illegal diag_style command");
      } else error->all(FLERR,"Illegal diag_style command");
    } else if (strcmp(arg[iarg],"delt") == 0) {
      iarg++;
      if (iarg < narg) {
	delta = atof(arg[iarg]);
	if (delta <= 0.0) error->all(FLERR,"Illegal diag_style command");
      } else error->all(FLERR,"Illegal diag_style command");
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal diag_style command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (scale <= 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat < 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"loglinfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal diag_style command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (scale <= 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat < 0) error->all(FLERR,"Illegal diag_style command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 2;
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      iarg++;
      if (iarg < narg) {
	delay = atof(arg[iarg]);
      } else error->all(FLERR,"Illegal diag_style command");
    } else break;
    iarg++;
  }

  iarg_child = iarg;

  if (stats_flag && logfreq) error->all(FLERR,"Illegal diag_style command");
  if (stats_flag && delta > 0.0) error->all(FLERR,"Illegal diag_style command");
  if (stats_flag && delay > 0.0) error->all(FLERR,"Illegal diag_style command");
}

/* ---------------------------------------------------------------------- */

Diag::~Diag()
{
  delete [] style;
}
