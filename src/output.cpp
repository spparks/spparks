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
#include "string.h"
#include "stdlib.h"
#include "output.h"
#include "style_dump.h"
#include "style_diag.h"
#include "app.h"
#include "dump.h"
#include "diag.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

#define MAXSTR 4096

/* ---------------------------------------------------------------------- */

Output::Output(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  stats_delta = 0.0;
  stats_logfreq = 0;
  stats_delay = 0.0;

  ndump = 0;
  max_dump = 0;
  dumplist = 0;

  ndiag = 0;
  diaglist = NULL;
}

/* ---------------------------------------------------------------------- */

Output::~Output()
{
  for (int i = 0; i < ndump; i++) delete dumplist[i];
  memory->sfree(dumplist);

  for (int i = 0; i < ndiag; i++) delete diaglist[i];
  memory->sfree(diaglist);
}

/* ---------------------------------------------------------------------- */

void Output::init(double time)
{
  for (int i = 0; i < ndump; i++) dumplist[i]->init();
  for (int i = 0; i < ndiag; i++) diaglist[i]->init();
}

/* ----------------------------------------------------------------------
   called before every run
   perform stats output
   set next output time for each kind of output
   return tnext = next time any output is needed
------------------------------------------------------------------------- */

double Output::setup(double time)
{
  // if dump has never occured, write initial snapshot
  // needs to happen in setup() in case propensity is output
  // b/c app not ready to compute propensities until setup_app() is called
  // set next time for each dump

  double dump_time = app->stoptime;
  for (int i = 0; i < ndump; i++) {
    if (dumplist[i]->idump == 0 && time >= dumplist[i]->delay)
      dumplist[i]->write(time);
    dumplist[i]->next_time = 
      next_time(time,dumplist[i]->logfreq,dumplist[i]->delta,
		dumplist[i]->nrepeat,dumplist[i]->scale,dumplist[i]->delay);
    dump_time = MIN(dump_time,dumplist[i]->next_time);
  }

  // if a diagnostic is drven by stats, compute the diagnostic
  // set next time for each diagnostic

  double diag_time = app->stoptime;
  for (int i = 0; i < ndiag; i++) {
    if  (diaglist[i]->stats_flag) diaglist[i]->compute();
    diaglist[i]->next_time = 
      next_time(time,diaglist[i]->logfreq,diaglist[i]->delta,
		diaglist[i]->nrepeat,diaglist[i]->scale,diaglist[i]->delay);
    diag_time = MIN(diag_time,diaglist[i]->next_time);
  }

  // perform stats output
  // set next time for stats

  stats_header();
  stats(0);
  stats_time = app->stoptime;
  if (stats_delta > 0.0)
    stats_time = next_time(time,stats_logfreq,stats_delta,
			   stats_nrepeat,stats_scale,stats_delay);

  // tnext = next output time for anything

  double tnext = app->stoptime;
  tnext = MIN(tnext,dump_time);
  tnext = MIN(tnext,diag_time);
  tnext = MIN(tnext,stats_time);
  return tnext;
}

/* ----------------------------------------------------------------------
   called only when some output is needed or when app is done
   set next output time for any output performed
   return tnext = next time any output is needed
------------------------------------------------------------------------- */

double Output::compute(double time, int done)
{
  // dump output

  double dump_time = app->stoptime;
  for (int i = 0; i < ndump; i++) {
    if (time >= dumplist[i]->next_time) {
      dumplist[i]->write(time);
      dumplist[i]->next_time = 
	next_time(time,dumplist[i]->logfreq,dumplist[i]->delta,
		  dumplist[i]->nrepeat,dumplist[i]->scale,dumplist[i]->delay);
      dump_time = MIN(dump_time,dumplist[i]->next_time);
    } else dump_time = MIN(dump_time,dumplist[i]->next_time);
  }

  // sflag = 1 if stats output needed

  int sflag = 0;
  if (time >= stats_time || done) sflag = 1;
  
  // diagnostic output, which may be driven by stats output
  
  double diag_time = app->stoptime;
  for (int i = 0; i < ndiag; i++) {
    if (diaglist[i]->stats_flag) {
      if (sflag) diaglist[i]->compute();
    } else if (time >= diaglist[i]->next_time) {
      diaglist[i]->compute();
      diaglist[i]->next_time = 
	next_time(time,diaglist[i]->logfreq,diaglist[i]->delta,
		  diaglist[i]->nrepeat,diaglist[i]->scale,diaglist[i]->delay);
      diag_time = MIN(diag_time,diaglist[i]->next_time);
    } else diag_time = MIN(diag_time,diaglist[i]->next_time);
  }
  
  // stats output, after diagnostics compute any needed quantities

  if (sflag) {
    stats(1);
    stats_time = app->stoptime;
    if (stats_delta)
      stats_time = next_time(time,stats_logfreq,stats_delta,
			     stats_nrepeat,stats_scale,stats_delay);
  }

  // tnext = next output time for anything

  double tnext = app->stoptime;
  tnext = MIN(tnext,dump_time);
  tnext = MIN(tnext,diag_time);
  tnext = MIN(tnext,stats_time);
  return tnext;
}

/* ---------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal stats command");
  stats_delta = atof(arg[0]);
  if (stats_delta < 0.0) error->all(FLERR,"Illegal stats command");
  if (stats_delta == 0.0 && narg > 1) error->all(FLERR,"Illegal stats command");

  stats_delay = 0.0;
  stats_logfreq = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal stats command");
      stats_delay = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal stats command");
      stats_nrepeat = atoi(arg[iarg+1]);
      stats_scale = atof(arg[iarg+2]);
      if (stats_scale <= 0) error->all(FLERR,"Illegal stats command");
      if (stats_nrepeat < 0) error->all(FLERR,"Illegal stats command");
      if (stats_nrepeat == 0) stats_logfreq = 0;
      else stats_logfreq = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"loglinfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal stats command");
      stats_nrepeat = atoi(arg[iarg+1]);
      stats_scale = atof(arg[iarg+2]);
      if (stats_nrepeat < 0) error->all(FLERR,"Illegal stats command");
      if (stats_nrepeat == 0) stats_logfreq = 0;
      else stats_logfreq = 2;
      iarg += 3;
    } else error->all(FLERR,"Illegal stats command");
  }
}

/* ---------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal dump command");

  for (int i = 0; i < ndump; i++)
    if (strcmp(dumplist[i]->id,arg[0]) == 0) 
      error->all(FLERR,"Reuse of dump ID");

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump++;
    dumplist = (Dump **)
      memory->srealloc(dumplist,max_dump*sizeof(Dump *),"output:dump");
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) dumplist[ndump] = new Class(spk,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  ndump++;
}

/* ----------------------------------------------------------------------
   force current snapshot to be written out
   does not change any attributes of next dump
------------------------------------------------------------------------- */

void Output::dump_one(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal dump_one command");

  int i;
  for (i = 0; i < ndump; i++)
    if (strcmp(dumplist[i]->id,arg[0]) == 0) break;
  if (i == ndump) 
    error->all(FLERR,"Could not find dump ID in dump_one command");

  if (dumplist[i]->idump == 0)
    error->all(FLERR,"Cannot use dump_one for first snapshot in dump file");

  dumplist[i]->write(app->time);
}

/* ---------------------------------------------------------------------- */

void Output::dump_modify(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal dump_modify command");

  int i;
  for (i = 0; i < ndump; i++)
    if (strcmp(dumplist[i]->id,arg[0]) == 0) break;
  if (i == ndump) 
    error->all(FLERR,"Could not find dump ID in dump_modify command");

  dumplist[i]->modify_params(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Output::undump(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal undump command");

  int i;
  for (i = 0; i < ndump; i++)
    if (strcmp(dumplist[i]->id,arg[0]) == 0) break;
  if (i == ndump) error->all(FLERR,"Could not find dump ID in undump command");

  delete dumplist[i];
  for (int j = i; j < ndump-1; j++) dumplist[j] = dumplist[j+1];
  ndump--;
}

/* ---------------------------------------------------------------------- */

void Output::add_diag(Diag *diag)
{
  ndiag++;
  diaglist = (Diag **) memory->srealloc(diaglist,ndiag*sizeof(Diag *),
					"output:diaglist");
  diaglist[ndiag-1] = diag;
}

/* ----------------------------------------------------------------------
   print stats, including contributions from app and diagnostics
------------------------------------------------------------------------- */

void Output::stats(int timeflag)
{
  char str[MAXSTR] = {'\0'};
  char *strpnt = str;

  app->stats(strpnt);
  strpnt += strlen(strpnt);

  if (timeflag) {
    sprintf(strpnt," %10.3g",timer->elapsed(TIME_LOOP));
    strpnt += strlen(strpnt);
  } else {
    sprintf(strpnt," %10.3g",0.0);
    strpnt += strlen(strpnt);
  }

  for (int i = 0; i < ndiag; i++)
    if (diaglist[i]->stats_flag) {
      diaglist[i]->stats(strpnt);
      strpnt += strlen(strpnt);
    }

  if (me == 0) {
    if (screen)
      fprintf(screen,"%s\n",str);
    if (logfile) {
      fprintf(logfile,"%s\n",str);
      fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   print stats header, including contributions from app and diagnostics
------------------------------------------------------------------------- */

void Output::stats_header()
{
  char str[MAXSTR] = {'\0'};
  char *strpnt = str;

  app->stats_header(strpnt);
  strpnt += strlen(strpnt);

  sprintf(strpnt," %10s","CPU");
  strpnt += strlen(strpnt);

  for (int i = 0; i < ndiag; i++)
    if (diaglist[i]->stats_flag) {
      diaglist[i]->stats_header(strpnt);
      strpnt += strlen(strpnt);
    }

  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",str);
    if (logfile) {
      fprintf(logfile,"%s\n",str);
      fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   calculate next time that output should be performed
   account for logarithmic frequency and delay
   return tnew = next time at which output should be done
------------------------------------------------------------------------- */

double Output::next_time(double tcurrent, int logfreq, double delta, 
			 int nrepeat, double scale, double delay)
{
  double tnew;

  if (logfreq == 0) {
    tnew = ceil(tcurrent/delta) * delta;
    if (tnew == tcurrent) tnew = tcurrent + delta;
  } else if (logfreq == 1) {
    while (tcurrent >= delta*scale) delta *= scale;
    double ktmp = pow(scale,1.0/nrepeat);
    tnew = delta;
    while (tcurrent >= tnew) tnew *= ktmp;
  } else if (logfreq == 2) {
    double start = delta;
    while (tcurrent >= start*scale) start *= scale;
    tnew = ceil(tcurrent/start) * start;
    if (tnew == tcurrent) tnew = tcurrent + start;
    if (static_cast<int> (tnew/start) > nrepeat) tnew = start*scale;
  }

  tnew = MAX(tnew,delay);
  return tnew;
}
