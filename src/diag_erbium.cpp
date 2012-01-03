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
#include "app_erbium.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{ZERO,ERBIUM,HYDROGEN,HELIUM,VACANCY};   // same as AppErbium
enum{ER,H,HE,VAC,EVENTS,ONE,TWO,THREE};

/* ---------------------------------------------------------------------- */

DiagErbium::DiagErbium(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"erbium") != 0)
    error->all(FLERR,"Diag_style erbium requires app_style erbium");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style erbium command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style erbium command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagErbium::~DiagErbium()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagErbium::init()
{
  apperbium = (AppErbium *) app;

  int none = apperbium->none;
  int ntwo = apperbium->ntwo;
  int nthree = apperbium->nthree;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"er") == 0) which[i] = ER;
    else if (strcmp(list[i],"h") == 0) which[i] = H;
    else if (strcmp(list[i],"he") == 0) which[i] = HE;
    else if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Invalid value setting in diag_style erbium");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Invalid value setting in diag_style erbium");
      index[i] = n - 1;
    } else if (list[i][0] == 't') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Invalid value setting in diag_style erbium");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style erbium");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
    if (which[i] == ER || which[i] == H || which[i] == HE || which[i] == VAC)
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagErbium::compute()
{
  int sites[5],ivalue;

  if (siteflag) {
    sites[ERBIUM] = sites[HYDROGEN] = sites[HELIUM] = sites[VACANCY] = 0;
    int *element = apperbium->element;
    int nlocal = apperbium->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == ER) ivalue = sites[ERBIUM];
    else if (which[i] == H) ivalue = sites[HYDROGEN];
    else if (which[i] == HE) ivalue = sites[HELIUM];
    else if (which[i] == VAC) ivalue = sites[VACANCY];
    else if (which[i] == EVENTS) ivalue = apperbium->nevents;
    else if (which[i] == ONE) ivalue = apperbium->scount[index[i]];
    else if (which[i] == TWO) ivalue = apperbium->dcount[index[i]];
    else if (which[i] == THREE) ivalue = apperbium->tcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErbium::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
