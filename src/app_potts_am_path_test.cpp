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

/* ---------------------------------------------------------------------- */

AppPottsAmPathTest::AppPottsAmPathTest(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg), xp(),yp()
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
   } else if (strcmp(command,"test_point") == 0) {
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      if (0==my_rank){
         printf("AppPottsAmPathTest::input_app \'test_point\'\n");
      }
      if (narg != 2) error->all(FLERR,"Illegal 'test_point' command; wrong num args; should have 2 arguments, x,y values.");
      double x=atof(arg[0]);
      double y=atof(arg[1]);
      xp.push_back(x);
      yp.push_back(y);
   } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAmPathTest::grow_app()
{
  spin = iarray[0];
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

   int num_test_points=xp.size();
	for(int i=0;i<num_test_points;i++){

		// Fake SPPARKS lattice site using input point defined by (xp,yp)
		double XYZ[]={xp[i],yp[i],0};

		// Lattice point location relative to 'pool' position
		Point xyz_r_p=compute_position_relative_to_pool(XYZ);

		//Temporary assignment of xo, xo is in the melt pool's reference frame!
		double xo[]={xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};

      // DEBUG
      //if (0==my_rank){
      //   std::cout << "Relative position: xo: " << std::setw(5) << xo[0] << "," << xo[1] << "," << xo[2] <<  std::endl;
      //}
      // END DEBUG
	}
}


/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
void AppPottsAmPathTest::site_event_rejection(int i, RandomPark *random) {
}
