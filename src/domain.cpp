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

#include "string.h"
#include "domain.h"
#include "app.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

#include "style_region.h"

using namespace SPPARKS_NS;

#define DELTA 1

/* ---------------------------------------------------------------------- */

Domain::Domain(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = 3;
  xperiodic = yperiodic = zperiodic = 1;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nx = ny = nz = 0;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;

  box_exist = 0;

  lattice = NULL;
  nregion = maxregion = 0;
  regions = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete lattice;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ----------------------------------------------------------------------
   setup global box
   assumes boxlo/hi are already set
------------------------------------------------------------------------- */

void Domain::set_box()
{
  if (boxxlo >= boxxhi || boxylo >= boxyhi || boxzlo >= boxzhi)
    error->one(FLERR,"Box bounds are invalid");

  xprd = boxxhi - boxxlo;
  yprd = boxyhi - boxylo;
  zprd = boxzhi - boxzlo;
}

/* ----------------------------------------------------------------------
   create a lattice
   delete it if style = none
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = new Lattice(spk,narg,arg);
  if (lattice->style == 0) {
    delete lattice;
    lattice = NULL;
  }
}

/* ----------------------------------------------------------------------
   create a new region 
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal region command");

  if (find_region(arg[0]) >= 0) error->all(FLERR,"Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTA;
    regions = (Region **) 
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");

#define REGION_CLASS
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    regions[nregion] = new Class(spk,narg,arg);
#include "style_region.h"
#undef REGION_CLASS

  else error->all(FLERR,"Invalid region style");

  nregion++;
}

/* ----------------------------------------------------------------------
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int Domain::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   boundary settings from the input script 
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  if (strcmp(arg[0],"n") == 0) xperiodic = 0;
  else if (strcmp(arg[0],"p") == 0) xperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[1],"n") == 0) yperiodic = 0;
  else if (strcmp(arg[1],"p") == 0) yperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");
  if (strcmp(arg[2],"n") == 0) zperiodic = 0;
  else if (strcmp(arg[2],"p") == 0) zperiodic = 1;
  else error->all(FLERR,"Illegal boundary command");

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) nonperiodic = 1;

  if (nonperiodic && app->appclass != App::LATTICE)
    error->all(FLERR,"Boundary command currently only supported by on-lattice apps");
}

/* ----------------------------------------------------------------------
   assign nprocs to 1d box as equal partitions
------------------------------------------------------------------------- */

void Domain::procs2domain_1d()
{
  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[1] != 1 || user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 1d simulation");
    procgrid[0] = user_procgrid[0];
  } else {
    procgrid[0] = nprocs;
  }

  procgrid[1] = procgrid[2] = 1;

  myloc[0] = me;
  myloc[1] = myloc[2] = 0;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo;
  subyhi = boxyhi;
  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void Domain::procs2domain_2d()
{
  int ipx,ipy;
  double boxx,boxy,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    if (user_procgrid[2] != 1)
      error->all(FLERR,"App style proc count is not valid for 2d simulation");
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];

  } else {

    // loop thru all possible factorizations of nprocs
    // surf = perimeter of a proc sub-domain

    double bestsurf = 2.0 * (xprd+yprd);
 
    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	ipy = nprocs/ipx;
	boxx = xprd/ipx;
	boxy = yprd/ipy;
	surf = boxx + boxy;
	if (surf < bestsurf) {
	  bestsurf = surf;
	  procgrid[0] = ipx;
	  procgrid[1] = ipy;
	}
      }
      ipx++;
    }
  }

  procgrid[2] = 1;

  myloc[0] = me % procgrid[0];
  myloc[1] = me/procgrid[0];
  myloc[2] = 0;

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1) 
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo;
  subzhi = boxzhi;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d box so as to minimize surface area per proc
------------------------------------------------------------------------- */

void Domain::procs2domain_3d()
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;

  if (user_procgrid[0] || user_procgrid[1] || user_procgrid[2]) {
    procgrid[0] = user_procgrid[0];
    procgrid[1] = user_procgrid[1];
    procgrid[2] = user_procgrid[2];

  } else {

    double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);
  
    // loop thru all possible factorizations of nprocs
    // surf = surface area of a proc sub-domain

    ipx = 1;
    while (ipx <= nprocs) {
      if (nprocs % ipx == 0) {
	nremain = nprocs/ipx;
	ipy = 1;
	while (ipy <= nremain) {
	  if (nremain % ipy == 0) {
	    ipz = nremain/ipy;
	    boxx = xprd/ipx;
	    boxy = yprd/ipy;
	    boxz = zprd/ipz;
	    surf = boxx*boxy + boxy*boxz + boxz*boxx;
	    if (surf < bestsurf) {
	      bestsurf = surf;
	      procgrid[0] = ipx;
	      procgrid[1] = ipy;
	      procgrid[2] = ipz;
	    }
	  }
	  ipy++;
	}
      }
      ipx++;
    }
  }

  myloc[0] = me % procgrid[0];
  myloc[1] = (me/procgrid[0]) % procgrid[1];
  myloc[2] = me / (procgrid[0]*procgrid[1]);

  subxlo = boxxlo + myloc[0] * xprd/procgrid[0];
  if (myloc[0] < procgrid[0]-1) 
    subxhi = boxxlo + (myloc[0]+1) * xprd/procgrid[0];
  else subxhi = boxxhi;

  subylo = boxylo + myloc[1] * yprd/procgrid[1];
  if (myloc[1] < procgrid[1]-1) 
    subyhi = boxylo + (myloc[1]+1) * yprd/procgrid[1];
  else subyhi = boxyhi;

  subzlo = boxzlo + myloc[2] * zprd/procgrid[2];
  if (myloc[2] < procgrid[2]-1) 
    subzhi = boxzlo + (myloc[2]+1) * zprd/procgrid[2];
  else subzhi = boxzhi;
}
