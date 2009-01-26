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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_pore.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED};

/* ---------------------------------------------------------------------- */

AppPore::AppPore(SPPARKS *spk, int narg, char **arg) : 
  AppDiffusionTable(spk,narg,arg)
{
  delevent = 1;
  delpropensity = 2;
  allow_rejection = 0;
  allow_masking = 0;

  // parse arguments

  if (narg < 7) error->all("Illegal app_style command");

  double xc = atof(arg[1]);
  double yc = atof(arg[2]);
  double zc = atof(arg[3]);
  double diameter = atof(arg[4]);
  double thickness = atof(arg[5]);
  int seed = atoi(arg[6]);
  random = new RandomPark(seed);

  options(narg-7,&arg[7]);

  // define lattice and partition it across processors
  // esites must be large enough for 2 sites and 1st/2nd nearest neighbors

  create_lattice();
  esites = new int[2 + 2*maxneigh + 2*maxneigh*maxneigh];
  echeck = NULL;

  // initialize my portion of lattice
  // each site = VACANT or OCCUPIED as defined by pore geometry

  if (infile) read_file();

  else {
    double x,y,z;
    int isite;
    for (int i = 0; i < nlocal; i++) {
      x = xyz[i][0];
      y = xyz[i][1];
      z = xyz[i][2];
      if (z > zc + 0.5*thickness || z < zc - 0.5*thickness) isite = VACANT;
      else isite = OCCUPIED;
      if (isite == OCCUPIED) {
	if ((x-xc)*(x-xc) + (y-yc)*(y-yc) < 0.25*diameter*diameter)
	  isite = VACANT;
      }
      lattice[i] = isite;
    }
  }
}
