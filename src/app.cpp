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
#include "stdlib.h"
#include "app.h"
#include "timer.h"
#include "finish.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

App::App(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  appclass = GENERAL;
  time = 0.0;
  first_run = 1;
}

/* ---------------------------------------------------------------------- */

App::~App()
{
  delete [] style;
}

/* ---------------------------------------------------------------------- */

void App::run(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal run command");

  stoptime = time + atof(arg[0]);
  if (stoptime < time) error->all("Illegal run command");

  // read optional args

  int uptoflag = 0;
  int preflag = 1;
  int postflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) error->all("Illegal run command");
      uptoflag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) preflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) preflag = 1;
      else error->all("Illegal run command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all("Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) postflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) postflag = 1;
      else error->all("Illegal run command");
      iarg += 2;
    } else error->all("Illegal run command");
  }

  // adjust stoptime if upto was specified

  if (uptoflag) stoptime -= time;
  if (uptoflag && stoptime < 0.0)
    error->all("Run upto value is before current time");

  // perform a single run via app's init(), setup(), and iterate()
  // if pre or 1st run, do app init
  // setup computes initial propensities
  // if post, do full Finish, else just print time

  if (preflag || first_run) {
    init();
    first_run = 0;
  }

  timer->init();
  setup();
  if (stoptime > time) iterate();

  Finish finish(spk,postflag);
}

/* ---------------------------------------------------------------------- */

void App::reset_time(double newtime)
{
  time = newtime;
}

/* ----------------------------------------------------------------------
   assign nprocs to 1d box as equal partitions
------------------------------------------------------------------------- */

void App::procs2domain_1d(int px_user, int py_user, int pz_user,
			  double xprd, double boxxlo, double boxxhi,
			  int &nx_procs,
			  double &subxlo, double &subylo, double &subzlo,
			  double &subxhi, double &subyhi, double &subzhi)
{
  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (px_user || py_user || pz_user) {
    if (py_user != 1 || pz_user != 1)
      error->all("App style proc count is not valid");
    nx_procs = px_user;
  } else {
    nx_procs = nprocs;
  }

  int iprocx = me;

  subxlo = boxxlo + iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = boxxlo + (iprocx+1) * xprd/nx_procs;
  else subxhi = boxxhi;

  subylo = -0.5;
  subyhi = 0.5;
  subzlo = -0.5;
  subzhi = 0.5;
}

/* ----------------------------------------------------------------------
   assign nprocs to 2d box so as to minimize perimeter per proc
------------------------------------------------------------------------- */

void App::procs2domain_2d(int px_user, int py_user, int pz_user,
			  double xprd, double yprd,
			  double boxxlo, double boxylo,
			  double boxxhi, double boxyhi,
			  int &nx_procs, int &ny_procs,
			  double &subxlo, double &subylo, double &subzlo,
			  double &subxhi, double &subyhi, double &subzhi)
{
  int ipx,ipy;
  double boxx,boxy,surf;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (px_user || py_user || pz_user) {
    if (pz_user != 1)
      error->all("App style proc count is not valid");
    nx_procs = px_user;
    ny_procs = py_user;

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
	  nx_procs = ipx;
	  ny_procs = ipy;
	}
      }
      ipx++;
    }
  }

  int iprocx = me/ny_procs;
  int iprocy = me % ny_procs;

  subxlo = boxxlo + iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = boxxlo + (iprocx+1) * xprd/nx_procs;
  else subxhi = boxxhi;

  subylo = boxylo + iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = boxylo + (iprocy+1) * yprd/ny_procs;
  else subyhi = boxyhi;

  subzlo = -0.5;
  subzhi = 0.5;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d box so as to minimize surface area per proc
------------------------------------------------------------------------- */

void App::procs2domain_3d(int px_user, int py_user, int pz_user,
			  double xprd, double yprd, double zprd,
			  double boxxlo, double boxylo, double boxzlo,
			  double boxxhi, double boxyhi, double boxzhi,
			  int &nx_procs, int &ny_procs, int &nz_procs,
			  double &subxlo, double &subylo, double &subzlo,
			  double &subxhi, double &subyhi, double &subzhi)
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (px_user || py_user || pz_user) {
    nx_procs = px_user;
    ny_procs = py_user;
    nz_procs = pz_user;

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
	      nx_procs = ipx;
	      ny_procs = ipy;
	      nz_procs = ipz;
	    }
	  }
	  ipy++;
	}
      }
      ipx++;
    }
  }

  int nyz_procs = ny_procs * nz_procs;

  int iprocx = (me/(nyz_procs)) % nx_procs;
  int iprocy = (me/nz_procs) % ny_procs;
  int iprocz = (me/1) % nz_procs;

  subxlo = boxxlo + iprocx * xprd/nx_procs;
  if (iprocx < nx_procs-1) subxhi = boxxlo + (iprocx+1) * xprd/nx_procs;
  else subxhi = boxxhi;

  subylo = boxylo + iprocy * yprd/ny_procs;
  if (iprocy < ny_procs-1) subyhi = boxylo + (iprocy+1) * yprd/ny_procs;
  else subyhi = boxyhi;

  subzlo = boxzlo + iprocz * zprd/nz_procs;
  if (iprocz < nz_procs-1) subzhi = boxzlo + (iprocz+1) * zprd/nz_procs;
  else subzhi = boxzhi;
}
