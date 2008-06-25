/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "spparks.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "input.h"
#include "app.h"
#include "solve.h"
#include "sweep.h"
#include "timer.h"
#include "output.h"

using namespace SPPARKS;

/* ----------------------------------------------------------------------
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

SPK::SPK(int narg, char **arg, MPI_Comm communicator)
{
  error = new Error(this);
  memory = new Memory(this);
  universe = new Universe(this,communicator);

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int iarg = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0) {
      if (iarg+1 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
	universe->add_world(arg[iarg]);
	iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0) {
      if (iarg+3 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg += 3;
    } else if (strcmp(arg[iarg],"-echo") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg += 2;
    } else error->universe_all("Invalid command-line argument");
  }

  // if procs was not a command-line switch, universe is one world w/ all procs

  if (universe->nworlds == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all("Processor partitions are inconsistent");

  // multiple-world universe must define input file

  if (universe->nworlds > 1 && inflag == 0)
    error->universe_all("Must use -in switch with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL) 
	error->universe_one("Cannot open universe screen file");
    }
    if (logflag == 0) {
      universe->ulogfile = fopen("log.spk","w");
      if (universe->ulogfile == NULL) 
	error->universe_one("Cannot open log.spk");
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL) 
	error->universe_one("Cannot open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // universe is single world
  // inherit settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->nworlds == 1) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;
    infile = NULL;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[inflag]);
	error->one(str);
      }
    }

    if (universe->me == 0) {
      if (screen) fprintf(screen,"SPPARKS (%s)\n",universe->version);
      if (logfile) fprintf(logfile,"SPPARKS (%s)\n",universe->version);
    }

  // universe is multiple worlds
  // split into separate communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    if (me == 0) {
      if (screenflag == 0) {
	char str[32];
	sprintf(str,"screen.%d",universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one("Cannot open screen file");
      } else if (strcmp(arg[screenflag],"none") == 0)
	screen = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one("Cannot open screen file");
      }
    } else screen = NULL;
    
    if (me == 0) {
      if (logflag == 0) {
	char str[32];
	sprintf(str,"log.spk.%d",universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one("Cannot open logfile");
      } else if (strcmp(arg[logflag],"none") == 0)
	logfile = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[logflag],universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one("Cannot open logfile");
      }
    } else logfile = NULL;
    
    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[inflag]);
	error->one(str);
      }
    } else infile = NULL;
    
    // screen and logfile messages for universe and world
    
    if (universe->me == 0) {
      if (universe->uscreen) {
	fprintf(universe->uscreen,"SPPARKS (%s)\n",universe->version);
	fprintf(universe->uscreen,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
      if (universe->ulogfile) {
	fprintf(universe->ulogfile,"SPPARKS (%s)\n",universe->version);
	fprintf(universe->ulogfile,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
    }
    
    if (me == 0) {
      if (screen) {
	fprintf(screen,"SPPARKS (%s)\n",universe->version);
	fprintf(screen,"Processor partition = %d\n",universe->iworld);
      }
      if (logfile) {
	fprintf(logfile,"SPPARKS (%s)\n",universe->version);
	fprintf(logfile,"Processor partition = %d\n",universe->iworld);
      }
    }
  }

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // allocate top-level classes

  create();
}

/* ----------------------------------------------------------------------
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

SPK::~SPK()
{
  destroy();

  if (universe->nworlds == 1) {
    if (logfile) fclose(logfile);
  } else {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
    if (universe->ulogfile) fclose(universe->ulogfile);
  }

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete input;
  delete universe;
  delete memory;
  delete error;
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in constructor
------------------------------------------------------------------------- */

void SPK::create()
{
  app = NULL;
  solve = NULL;
  sweep = NULL;
  timer = new Timer(this);
  output = new Output(this);
}

/* ----------------------------------------------------------------------
   init top-level classes
------------------------------------------------------------------------- */

void SPK::init()
{
  //app->init();
  //solve->init();
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void SPK::destroy()
{
  delete app;
  delete solve;
  delete sweep;
  delete timer;
  delete output;
}
