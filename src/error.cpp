/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "error.h"
#include "universe.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Error::Error(SPPARKS *spk) : Pointers(spk) {}

/* ----------------------------------------------------------------------
   called by all procs in universe
   close all output, screen, and log files in world and universe
   no abort, so insure all procs in universe call, else will hang
------------------------------------------------------------------------- */

void Error::universe_all(const char *file, int line, const char *str)
{
  MPI_Barrier(universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,
				   "ERROR: %s (%s:%d)\n",str,file,line);
    if (universe->ulogfile) fprintf(universe->ulogfile,
				    "ERROR: %s (%s:%d)\n",str,file,line);
  }

  //if (output) delete output;
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   forces abort of entire universe if any proc in universe calls
------------------------------------------------------------------------- */

void Error::universe_one(const char *file, int line, const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
	    universe->me,str,file,line);
  MPI_Abort(universe->uworld,1);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   insure all procs in world call, else will hang
   force MPI_Abort if running in multi-partition mode
------------------------------------------------------------------------- */

void Error::all(const char *file, int line, const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if (screen) fprintf(screen,"ERROR: %s (%s:%d)\n",str,file,line);
    if (logfile) fprintf(logfile,"ERROR: %s (%s:%d)\n",str,file,line);
  }

  //if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-NULL on this proc
   always write to universe screen 
   forces abort of entire world (and universe) if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const char *file, int line, const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);
  if (screen) fprintf(screen,"ERROR on proc %d: %s (%s:%d)\n",
		      me,str,file,line);
  if (universe->nworlds > 1)
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
	    universe->me,str,file,line);
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file 
------------------------------------------------------------------------- */

void Error::warning(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"WARNING: %s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"WARNING: %s (%s:%d)\n",
				  str,file,line);
}

/* ----------------------------------------------------------------------
   called by one proc in world, typically proc 0
   write message to screen and logfile (if logflag is set)
------------------------------------------------------------------------- */

void Error::message(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"%s (%s:%d)\n",str,file,line);
  if (logflag && logfile) fprintf(logfile,"%s (%s:%d)\n",str,file,line);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   no abort, so insure all procs in world call, else will hang
------------------------------------------------------------------------- */

void Error::done()
{
  MPI_Barrier(world);

  //if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  MPI_Finalize();
  exit(1);
}
