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
  stats_delta = 0.0;
  dump_delta = 0.0;
  ndiags = 0;
  diaglist = NULL;
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

  stats_time = time + stats_delta;
  dump_time = time + dump_delta;

  // print dump file header and 1st snapshot

  if (dump_delta > 0.0) {
    app->dump_header();
    app->dump();
  }

  // print stats header and initial stats
  
  app->stats_header();
  app->stats();

  // Initialize all diagnostics

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->init(time);
  }

}


/* ---------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
  app->set_stats(narg, arg);
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
  if ((stats_delta > 0.0 && time >= stats_time) || done) {
    app->stats();
    stats_time += stats_delta;
    timer->stamp(TIME_OUTPUT);
  }
  
  if ((dump_delta > 0.0 && time >= dump_time) || done) {
    if (dump_delta > 0.0) app->dump();
    dump_time += dump_delta;
    timer->stamp(TIME_OUTPUT);
  }

  // Perform all diagnostics

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->compute(time,done);
  }


}

/* ---------------------------------------------------------------------- */

void Output::add_diag(Diag *diag)
{
  ndiags++;
  diaglist = (Diag **) memory->srealloc(diaglist,ndiags*sizeof(Diag*),"output:diaglist");
  diaglist[ndiags-1] = diag;
}

