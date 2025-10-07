/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "spktype.h"
#include "mpi.h"
#include "string.h"
#include "spparks.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "input.h"
#include "app.h"
#include "solve.h"
#include "domain.h"
#include "potential.h"
#include "random_mars.h"
#include "timer.h"
#include "output.h"

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

SPPARKS::SPPARKS(int narg, char **arg, MPI_Comm communicator)
{
  error = new Error(this);
  memory = new Memory(this);
  universe = new Universe(this,communicator);

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int helpflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0) {
      if (iarg+1 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
	universe->add_world(arg[iarg]);
	iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0) {
      if (iarg+3 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 3;
    } else if (strcmp(arg[iarg],"-echo") == 0) {
      if (iarg+2 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;
    } else if (strcmp(arg[iarg],"-help") == 0 || 
	       strcmp(arg[iarg],"-h") == 0) {
      if (iarg+1 > narg) 
	error->universe_all(FLERR,"Invalid command-line argument");
      helpflag = 1;
      iarg += 1;
    } else error->universe_all(FLERR,"Invalid command-line argument");
  }

  // if procs was not a command-line switch, universe is one world w/ all procs

  if (universe->nworlds == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all(FLERR,"Processor partitions are inconsistent");

  // multiple-world universe must define input file

  if (universe->nworlds > 1 && inflag == 0)
    error->universe_all(FLERR,"Must use -in switch with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL) 
	error->universe_one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      universe->ulogfile = fopen("log.spparks","w");
      if (universe->ulogfile == NULL) 
	error->universe_one(FLERR,"Cannot open log.spparks");
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL) 
	error->universe_one(FLERR,"Cannot open universe log file");
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
	error->one(FLERR,str);
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
	if (screen == NULL) error->one(FLERR,"Cannot open screen file");
      } else if (strcmp(arg[screenflag],"none") == 0)
	screen = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one(FLERR,"Cannot open screen file");
      }
    } else screen = NULL;
    
    if (me == 0) {
      if (logflag == 0) {
	char str[32];
	sprintf(str,"log.spparks.%d",universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      } else if (strcmp(arg[logflag],"none") == 0)
	logfile = NULL;
      else {
	char str[128];
	sprintf(str,"%s.%d",arg[logflag],universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one(FLERR,"Cannot open logfile");
      }
    } else logfile = NULL;
    
    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) {
	char str[128];
	sprintf(str,"Cannot open input script %s",arg[inflag]);
	error->one(FLERR,str);
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

  // check datatype settings in spparks.h

  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in spktype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in spktype.h is invalid");
  if (sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in spktype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_SPK_TAGINT,&mpisize);
  if (mpisize != sizeof(tagint))
      error->all(FLERR,
		 "MPI_SPK_TAGINT and tagint in spktype.h are not compatible");
  MPI_Type_size(MPI_SPK_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
      error->all(FLERR,
		 "MPI_SPK_BIGINT and bigint in spktype.h are not compatible");

  // allocate input class now that MPI is fully setup

  input = new Input(this,narg,arg);

  // allocate top-level classes

  create();

  // if helpflag set, print help and exit

  if (helpflag) {
    if (universe->me == 0) print_styles();
    error->done();
  }
}

/* ----------------------------------------------------------------------
   shutdown SPPARKS
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

SPPARKS::~SPPARKS()
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

void SPPARKS::create()
{
  app = NULL;
  solve = NULL;
  domain = new Domain(this);
  potential = new Potential(this);
  ranmaster = new RanMars(this);
  output = new Output(this);
  timer = new Timer(this);
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void SPPARKS::destroy()
{
  delete app;
  delete solve;
  delete domain;
  delete potential;
  delete ranmaster;
  delete output;
  delete timer;
}

/* ----------------------------------------------------------------------
   for each style, print name of all child classes build into executable
------------------------------------------------------------------------- */

void SPPARKS::print_styles()
{
  printf("\nList of style options included in this executable:\n\n");

  printf("App styles:");
#define APP_CLASS
#define AppStyle(key,Class) printf(" %s",#key);
#include "style_app.h"
#undef APP_CLASS
  printf("\n\n");

  printf("Solve styles:");
#define SOLVE_CLASS
#define SolveStyle(key,Class) printf(" %s",#key);
#include "style_solve.h"
#undef SOLVE_CLASS
  printf("\n\n");

  printf("Diag styles:");
#define DIAG_CLASS
#define DiagStyle(key,Class) printf(" %s",#key);
#include "style_diag.h"
#undef DIAG_CLASS
  printf("\n\n");

  printf("Dump styles:");
#define DUMP_CLASS
#define DumpStyle(key,Class) printf(" %s",#key);
#include "style_dump.h"
#undef DUMP_CLASS
  printf("\n\n");

  printf("Pair styles:");
#define PAIR_CLASS
#define PairStyle(key,Class) printf(" %s",#key);
#include "style_pair.h"
#undef PAIR_CLASS
  printf("\n\n");

  printf("Region styles:");
#define REGION_CLASS
#define RegionStyle(key,Class) printf(" %s",#key);
#include "style_region.h"
#undef REGION_CLASS
  printf("\n\n");

  printf("Command styles (add-on input script commands):");
#define COMMAND_CLASS
#define CommandStyle(key,Class) printf(" %s",#key);
#include "style_command.h"
#undef COMMAND_CLASS
  printf("\n");
}
