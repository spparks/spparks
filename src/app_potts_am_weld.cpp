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

#include <iostream>
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
#include "am_teardrop.h"
#include "app_potts_am_weld.h"

using namespace SPPARKS_NS;

namespace RASTER {

namespace pool_shape {

namespace AM_TEARDROP {

namespace DEMO {

   enum InchesPerMinute {I=20,II=25,III=30};

   vector< vector<double> > get_demo_control_points(InchesPerMinute ipm) {
       int n=5;
       vector< vector<double> > cp(n);
       for(int i=0;i<n;i++)
          cp[i]=vector<double>(2);
       double a[]={0.0,0.0}, b[]={0.0,0.0}, c[]={0.0,0.0}, d[]={0.0,0.0}, e[]={0.0,0.0,};
       switch(ipm){
          case I:
             a[0]=-6.35;  a[1]=0.00;
             b[0]=a[0];   b[1]=1.50;
             c[0]=-3.25;  c[1]=3.00;
             d[0]=0.00;   d[1]=2.75;
             e[0]=d[0];   e[1]=0.00;
             break;
          case II:
             a[0]=-7.35;  a[1]=0.00;
             b[0]=a[0];   b[1]=1.00;
             c[0]=-3.75;  c[1]=2.60;
             d[0]=0.00;   d[1]=2.80;
             e[0]=d[0];   e[1]=0.00;
             break;
          case III:
             a[0]=-7.35;  a[1]=0.00;
             b[0]=a[0];   b[1]=0.70;
             c[0]=-3.75;  c[1]=2.00;
             d[0]=0.00;   d[1]=2.475;
             e[0]=d[0];   e[1]=0.00;
             break;
       }
       cp[0][0]= a[0];  cp[0][1]= a[1];
       cp[1][0]= b[0];  cp[1][1]= b[1];
       cp[2][0]= c[0];  cp[2][1]= c[1];
       cp[3][0]= d[0];  cp[3][1]= d[1];
       cp[4][0]= e[0];  cp[4][1]= e[1];
       return cp;
   }

} // DEMO

} // AM_TEARDROP


} // pool_shape

} // RASTER

using RASTER::pool_shape::AM_TEARDROP::DEMO::InchesPerMinute;
using RASTER::pool_shape::AM_TEARDROP::DEMO::get_demo_control_points;
using RASTER::pool_shape::ShapeType;

/* ---------------------------------------------------------------------- */

AppPottsAmWeld::AppPottsAmWeld(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg), alpha(0.5), beta(0.75), haz(-1.0), distance(nullptr),
   random_park(std::atof(arg[2])), simulation_time(0.0), shape_type(ShapeType::undefined),
   width(-1.0), length(-1.0), teardrop_control_points()
{

   // only error check for this class, not derived classes
   if (std::strcmp(arg[0],"potts/am/weld") != 0 || narg != 5 )
      error->all(FLERR,"Illegal app_style in 'potts/am/weld' command");

   // Flag which forces 'callback' to this app each step time 'time' is updated;
   // See 'allow_app_update' in app_lattice.h
   allow_app_update=1;
   // Adding 'distance' array; number of 'double' values per site
   // app.h
   ndouble = 1;

   // app_potts.cpp
   nspins = atoi(arg[1]);

   // app_potts_am_weld model parameters
   alpha =     std::atof(arg[2]); // relative size of pool shape at bottom compared to top
   beta =      std::atof(arg[3]); // defines curvature of pool shape through thickness
   haz =       std::atof(arg[4]); // distance which defines 'heat affected zone'

   // add the double array
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAmWeld::input_app(char *command, int narg, char **arg)
{   
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else if(strcmp(command,"weld_shape_teardrop")==0) {
      shape_type=RASTER::pool_shape::ShapeType::teardrop;
      if(strcmp(arg[0],"width")==0){
         width =  std::atof(arg[1]); // weld pool width (y-coordinate axis)
      } else error->all(FLERR,"Unrecognized 'weld_shape_teardrop' command.  Expected 'width'.");
      if(strcmp(arg[2],"case")==0){
         if(strcmp("I",arg[3])==0) teardrop_control_points=get_demo_control_points(InchesPerMinute::I);
         else if(strcmp("II",arg[3])==0) teardrop_control_points=get_demo_control_points(InchesPerMinute::II);
         else if(strcmp("III",arg[3])==0) teardrop_control_points=get_demo_control_points(InchesPerMinute::III);
         else if(strcmp("control_points",arg[3])==0) error->all(FLERR,"NOT IMPLEMENTED 'weld_shape_teardrop control_points' command.");
         else error->all(FLERR,"Unrecognized 'weld_shape_teardrop case' command");
      } else error->all(FLERR,"Unrecognized 'weld_shape_teardrop' command");
   } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAmWeld::grow_app()
{
  spin = iarray[0];
  distance = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAmWeld::init_app()
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

void AppPottsAmWeld::app_update(double dt)
{
   // Move pool
   bool moved=app_update_am(dt);
   if(!moved)
      return;

   // HACK
   // TODO: FIX ME
   const Domain *domain=this->domain;
   double lz=std::abs(domain->zprd);
   double thickness=lz;
   // std::cout << "Plate thickness: " << std::setw(5) << thickness << std::endl;
   // END HACK

   RASTER::pool_shape::PoolShape *shape=nullptr;
   switch(shape_type) {
      case RASTER::pool_shape::ShapeType::teardrop:
         shape_type=RASTER::pool_shape::ShapeType::teardrop;
         shape=new RASTER::pool_shape::AM_TEARDROP::Teardrop(thickness,width,alpha,beta,teardrop_control_points);
         break;
   }

	//Go through all the local sites and calculate the distance.
	for(int i=0;i<nlocal;i++){

		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};
		// Lattice point location relative to 'pool' position
		Point xyz_r_p=compute_position_relative_to_pool(XYZ);

		//Temporary assignment of xo, xo is in the melt pool's reference frame!
		double xo[]={xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};

      // Compute distance for this site
      // This function returns d=-1 for points within pool;
      // Otherwise it returns the distance associated with the closest point
      double d=shape->distance(xo);
      distance[i]=d;

      /* -------------------------------------------------------------------
      For sites within pool, assign random spins so that when they exit the
      pool, they will be singleton spins and generally not part of a grain.
      ---------------------------------------------------------------------- */
      if(d<0){
         // Random number between [1,nspins] inclusive
         int ran = random_park.irandom(nspins);
         spin[i]=ran;
      }

	}
}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double AppPottsAmWeld::compute_mobility(int site)  {
   double d=distance[site];
   //   mobility = 0 @ d<0
   //   mobility = 1 @ d=0 pool/solid interface
   //   mobility = 0 @ d>=haz distance defining heat effected zone
   //   Otherwise, mobility decreases linearly for 0<=d<haz
   double mobility=1-d/haz;
   if (d<0 || d>=haz)
      mobility=0.0;
   return mobility;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsAmWeld::site_event_rejection(int i, RandomPark *random) {
   int oldstate = spin[i];
   double einitial = site_energy(i);
   double mobility=compute_mobility(i);
   int nevent = 0;

   if(mobility>0.0){
      /*
       * Finding unique neighbor spins
       */
      std::set<int> my_unique;
      for (int j = 0; j < numneigh[i]; j++) {
         int value = spin[neighbor[i][j]];
         double j_distance = distance[neighbor[i][j]];
         // Don't flip to 'spins' on interior since they are not 'grains'
         // distance <= 0 are 'pinned' sites on the interior of weld pool;
         if (value == spin[i] || j_distance < 0.0) continue;
         my_unique.insert(value);
      }
      // how many unique spins
      nevent=my_unique.size();

      // scale random number according number of unique neighbor spins
      if (nevent == 0) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;

      {
         /*
          * Store unique spins
          */
         int c=0;
          for(std::set<int>::iterator i=my_unique.begin();i!=my_unique.end();i++){
            unique[c++]=*i;
          }
      }
      spin[i] = unique[iran];
      double efinal = site_energy(i);

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > mobility){
               spin[i] = oldstate;
         }
      } else if (0.0==temperature) {
         spin[i] = oldstate;
      } else if (random->uniform() > mobility * exp((einitial-efinal)*t_inverse)) {
         spin[i] = oldstate;
      }
  }

  if (spin[i] != oldstate) naccept++;
}
