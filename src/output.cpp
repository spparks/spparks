/* -----------------------q-----------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

Output::Output(SPK *spk) : SysPtr(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  stats_delta = 0.0;
  dump_delta = 0.0;
  ndiags = 0;
  diaglist = NULL;
  stats_ilogfreq = 0;
}

/* ---------------------------------------------------------------------- */

Output::~Output()
{
  // Free memory allocated to diagnostics.
  for (int i = 0; i < ndiags; i++) {
    delete diaglist[i];
  }
  memory->sfree(diaglist);
}

/* ---------------------------------------------------------------------- */

void Output::init(double time)
{
  // setup future stat and dump calls

  dump_time = time + dump_delta;

  if (stats_ilogfreq == 0) {
    stats_time = time + stats_delta;
  } else if (stats_ilogfreq == 1) {
    stats_time = time + stats_delta;
    stats_t0 = time;
    stats_irepeat = 0;
  }

  // print dump file header and 1st snapshot

  if (dump_delta > 0.0) {
    app->dump_header();
    app->dump();
  }

  // Initialize all diagnostics

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->init(time);
  }

  // print stats header and initial stats
  
  stats_header();
  stats();

}


/* ---------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      stats_ilogfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	stats_nrepeat = atoi(arg[iarg]);
	iarg++;
	stats_scale = atof(arg[iarg]);
      } else {
	error->all("Illegal stats command");
      }
    }
    iarg++;
  }

  app->set_stats(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Output::set_dump(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  app->set_dump(narg, arg);
}


/* ---------------------------------------------------------------------- */

void Output::compute(double time, int done)
{
  if ((dump_delta > 0.0 && time >= dump_time) || done) {
    if (dump_delta > 0.0) app->dump();
    dump_time += dump_delta;
  }

  // Perform all diagnostics

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->compute(time,done);
  }

  if ((stats_delta > 0.0 && time >= stats_time) || done) {
    if (stats_ilogfreq == 0) {
      stats_time += stats_delta;
    } else if (stats_ilogfreq == 1) {
      stats_time += stats_delta;
      stats_irepeat++;
      if (stats_irepeat == stats_nrepeat) {
	stats_delta *= stats_scale;
	stats_time = stats_t0+stats_delta;
	stats_irepeat = 0;
      }
    }
    stats();
  }
  
}

/* ---------------------------------------------------------------------- */

void Output::add_diag(Diag *diag)
{
  ndiags++;
  diaglist = (Diag **) memory->srealloc(diaglist,ndiags*sizeof(Diag *),"output:diaglist");
  diaglist[ndiags-1] = diag;
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void Output::stats()
{
  char str[1024] = {'\0'};
  char *strpnt = str;

  app->stats(strpnt);
  strpnt += strlen(strpnt);

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->stats(strpnt);
    strpnt += strlen(strpnt);
  }
  if (me == 0) {
    if (screen)
      fprintf(screen,"%s\n",str);
    if (logfile)
      fprintf(logfile,"%s\n",str);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void Output::stats_header()
{
  char str[1024] = {'\0'};
  char *strpnt = str;

  app->stats_header(strpnt);
  strpnt += strlen(strpnt);

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->stats_header(strpnt);
    strpnt += strlen(strpnt);
  }

  if (me == 0) {
    if (screen)
      fprintf(screen,"%s\n",str);
    if (logfile)
      fprintf(logfile,"%s\n",str);
  }
}
