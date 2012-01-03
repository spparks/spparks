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
#include "domain.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
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

  ninteger = ndouble = 0;
  id = NULL;
  xyz = NULL;
  iarray = NULL;
  darray = NULL;

  sites_exist = 0;
}

/* ---------------------------------------------------------------------- */

App::~App()
{
  delete [] style;

  memory->destroy(id);
  memory->destroy(xyz);
  for (int i = 0; i < ninteger; i++) memory->destroy(iarray[i]);
  for (int i = 0; i < ndouble; i++) memory->destroy(darray[i]);
  delete [] iarray;
  delete [] darray;
}

/* ---------------------------------------------------------------------- */

void App::create_arrays()
{
  if (ninteger) iarray = new int*[ninteger];
  for (int i = 0; i < ninteger; i++) iarray[i] = NULL;
  if (ndouble) darray = new double*[ndouble];
  for (int i = 0; i < ndouble; i++) darray[i] = NULL;
}

/* ---------------------------------------------------------------------- */

void App::recreate_arrays()
{
  delete [] iarray;
  delete [] darray;
  create_arrays();
}

/* ---------------------------------------------------------------------- */

void App::run(int narg, char **arg)
{
  if (appclass != GENERAL && domain->box_exist == 0)
    error->all(FLERR,"Cannot run application until simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal run command");

  stoptime = time + atof(arg[0]);
  if (stoptime < time) error->all(FLERR,"Illegal run command");

  // read optional args

  int uptoflag = 0;
  int preflag = 1;
  int postflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal run command");
      uptoflag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) preflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) preflag = 1;
      else error->all(FLERR,"Illegal run command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal run command");
      if (strcmp(arg[iarg+1],"no") == 0) postflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) postflag = 1;
      else error->all(FLERR,"Illegal run command");
      iarg += 2;
    } else error->all(FLERR,"Illegal run command");
  }

  // adjust stoptime if upto was specified

  if (uptoflag) stoptime -= time;
  if (uptoflag && stoptime < 0.0)
    error->all(FLERR,"Run upto value is before current time");

  // perform a single run via app's init(), setup(), and iterate()
  // if pre or 1st run, do app init
  // setup computes initial propensities
  // if post, do full Finish, else just print time

  if (preflag || first_run) {
    init();
    first_run = 0;
  }

  if (domain->me == 0) {
    if (screen) fprintf(screen,"Setting up run ...\n");
    if (logfile) fprintf(logfile,"Setting up run ...\n");
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
   return a pointer to a named internal variable
   if don't recognize name, pass it along to lower-level app
   names iarrayN and darrayN mean entry N from 1 to ninteger or ndouble
 ------------------------------------------------------------------------- */

void *App::extract(char *name)
{
  if (strcmp(name,"dimension") == 0) return (void *) &domain->dimension;
  if (strcmp(name,"boxxlo") == 0) return (void *) &domain->boxxlo;
  if (strcmp(name,"boxxhi") == 0) return (void *) &domain->boxxhi;
  if (strcmp(name,"boxylo") == 0) return (void *) &domain->boxylo;
  if (strcmp(name,"boxyhi") == 0) return (void *) &domain->boxyhi;
  if (strcmp(name,"boxzlo") == 0) return (void *) &domain->boxzlo;
  if (strcmp(name,"boxzhi") == 0) return (void *) &domain->boxzhi;

  if (strcmp(name,"nglobal") == 0) return (void *) &nglobal;
  if (strcmp(name,"nlocal") == 0) return (void *) &nlocal;

  if (strcmp(name,"id") == 0) return (void *) id;
  if (strcmp(name,"xyz") == 0) return (void *) xyz;

  if (strcmp(name,"site") == 0) {
    if (ninteger == 0) return NULL;
    return (void *) iarray[0];
  }

  if (strstr(name,"iarray") == name) {
    int n = atoi(&name[6]);
    if (n < 1 || n > ninteger) return NULL;
    return (void *) iarray[n-1];
  }
  if (strstr(name,"darray") == name) {
    int n = atoi(&name[6]);
    if (n < 1 || n > ndouble) return NULL;
    return (void *) darray[n-1];
  }

  return extract_app(name);
}

/* ----------------------------------------------------------------------
   return max ID
   may not be nglobal if site IDs are not contiguous
 ------------------------------------------------------------------------- */

tagint App::max_site_ID()
{
  tagint max = -1;
  for (int i = 0; i < nlocal; i++) max = MAX(max,id[i]);
  tagint all;
  MPI_Allreduce(&max,&all,1,MPI_SPK_TAGINT,MPI_MAX,world);
  return all;
}

/* ----------------------------------------------------------------------
   check if site IDs are contiguous from 1 to N
 ------------------------------------------------------------------------- */

int App::contiguous_sites()
{
  tagint min = nglobal+1;
  tagint max = -1;

  for (int i = 0; i < nlocal; i++) {
    min = MIN(min,id[i]);
    max = MAX(max,id[i]);
  }

  tagint all;
  MPI_Allreduce(&min,&all,1,MPI_SPK_TAGINT,MPI_MIN,world);
  if (all != 1) return 0;
  MPI_Allreduce(&max,&all,1,MPI_SPK_TAGINT,MPI_MAX,world);
  if (all != nglobal) return 0;
  return 1;
}

