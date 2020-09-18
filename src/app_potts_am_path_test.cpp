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

#include <stdio.h>
#include <cstdlib>
#include <cstddef>
#include <string>
#include <cstring>
#include <cmath>
#include <set>
#include <vector>
#include "string.h"
#include "math.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "error.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include "app_potts_am_path_test.h"

using namespace SPPARKS_NS;
using RASTER::DIR;

/* ---------------------------------------------------------------------- */

AppPottsAmPathTest::AppPottsAmPathTest(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg)
{
   // only error check for this class, not derived classes
   if (std::strcmp(arg[0],"potts/am/path/test") != 0 || narg != 2 )
      error->all(FLERR,"Illegal app_style in 'potts/am/path/test' command");

   // Flag which forces 'callback' to this app each step time 'time' is updated;
   // See 'allow_app_update' in app_lattice.h
   allow_app_update=1;

   // app_potts.cpp
   nspins = atoi(arg[1]);
//   int my_rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//   if (0==my_rank){
//      printf("%s\n", "AppPottsAmPathTest::AppPottsAmPathTest() ");
//   }
}

/* ---------------------------------------------------------------------- */

AppPottsAmPathTest::~AppPottsAmPathTest() { 
//   int my_rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//   if (0==my_rank){
//      printf("%s\n", "AppPottsAmPathTest::~AppPottsAmPathTest() DESTRUCTOR called. ");
//   }
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAmPathTest::input_app(char *command, int narg, char **arg)
{
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else if (strcmp(command,"pathgen") == 0) {
      if (narg != 12) error->all(FLERR,"Illegal 'pathgen' command; wrong num args; should have 8 arguments.");
      if(strcmp(arg[0],"outfile")==0){
         char *filename=arg[1];
         path_filename=std::string(filename);
      } else error->all(FLERR,"Unrecognized command: Expected 'outfile'.");
      if(strcmp(arg[2],"num_layers")==0){
         num_layers=atoi(arg[3]);
      }else error->all(FLERR,"Unrecognized command: Expected 'num_layers'.");
      if(strcmp(arg[4],"zstart")==0){
         zstart=atoi(arg[5]);
      }else error->all(FLERR,"Unrecognized command: Expected 'zstart'.");
      if(strcmp(arg[6],"width_haz")==0){
         width_haz=atoi(arg[7]);
      }else error->all(FLERR,"Unrecognized command: Expected 'width_haz'.");
      if(strcmp(arg[8],"melt_depth")==0){
         melt_depth=atoi(arg[9]);
      }else error->all(FLERR,"Unrecognized command: Expected 'melt_depth'.");
      if(strcmp(arg[10],"depth_haz")==0){
         depth_haz=atoi(arg[11]);
      }else error->all(FLERR,"Unrecognized command: Expected 'depth_haz'.");
      // 'pathgen' command must come after all 'am' commands in script
      print_paths(path_filename,num_layers,melt_depth,width_haz,depth_haz);
   } else error->all(FLERR,"Unrecognized command");
   // print path
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAmPathTest::init_app()
{
   // Run base class init_app
   init_app_am();
   // Compute distance function based upon initial pool position
   app_update(0.0);

}

/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */

void AppPottsAmPathTest::app_update(double dt)
{
   // Move pool
   bool moved=app_update_am(dt);
   if(!moved)
      return;
}
