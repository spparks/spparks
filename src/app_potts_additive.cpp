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

/* ----------------------------------------------------------------------
   Contributing author: Theron Rodgers and John Mitchell (Sandia)
------------------------------------------------------------------------- */

#include "string.h"
#include "math.h"
#include "app_potts_additive.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "am_ellipsoid.h"
#include "error.h"

#include <iostream>

using namespace SPPARKS_NS;
using RASTER::point;
using RASTER::DIR;

/* ---------------------------------------------------------------------- */

AppPottsAdditive::AppPottsAdditive(SPPARKS *spk, int narg, char **arg) :
  AppPotts(spk,narg,arg), passes(), active_layer() {

   // only error check for this class, not derived classes
   if (strcmp(arg[0],"additive") == 0 && narg != 11 )
    error->all(FLERR,"Illegal app_style command");
    
   nspins = atoi(arg[1]); //Number of spins
   spot_width = atof(arg[2]); //Width of the melt pool
   melt_tail_length = atof(arg[3]); //Length of tail from meltpool midpoint
   melt_depth = atof(arg[4]); //How many lattice sites deep the melt pool is
   cap_height = atof(arg[5]); //Height of the cap leading the meltpool
   HAZ = atof(arg[6]); //Size of the HAZ surrounding the melt pool (must be larger than spot_width)
   tail_HAZ = atof(arg[7]); //Length of hot zone behind meltpool (must be larger than melt_tail_length)
   depth_HAZ = atof(arg[8]); //Depth of the hot zone underneath the meltpool (must be larger than melt_depth)
   cap_HAZ = atof(arg[9]); //Size of HAZ infront of the melt pool (must be larger than cap_height)
   exp_factor = atof(arg[10]); //Exponential parameter for mobility decay in haz M(d) = exp(-exp_factor * d)

   //Define the layer object, this might work better in init_app
   ndouble = 1;
   allow_app_update = 1;
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAdditive::input_app(char *command, int narg, char **arg) 
{
   if (strcmp(command,"am_pass") == 0) {
      if (narg != 7) error->all(FLERR,"Illegal pass command.");
      char* dir=nullptr;
      double distance=-1.0;
      double speed=-1.0;
      int id=std::atoi(arg[0]);
      if(strcmp(arg[1],"dir")==0){
         dir=arg[2];
      } else {error->all(FLERR,"Illegal pass command. Expected 'dir.'");}
      if(strcmp(arg[3],"distance")==0){
         distance=std::atof(arg[4]);
      } else {error->all(FLERR,"Illegal pass command. Expected 'distance.'");}
      if(strcmp(arg[5],"speed")==0){
         speed=std::atof(arg[6]);
      } else {error->all(FLERR,"Illegal pass command. Expected 'speed.'");}
      DIR d;
      if(strcmp(dir,"X")==0) d=DIR::X;
      else if(strcmp(dir,"Y")==0) d=DIR::Y;
      else {error->all(FLERR,"Illegal pass 'dir' command. Expected 'X|Y.'");}
      passes[id]=Pass(d,distance,speed);

   } else if (strcmp(command,"am_transverse_pass") == 0) {
      if (narg != 5) error->all(FLERR,"Illegal transverse_pass command.");
      double distance=-1.0;
      double increment=-1.0;
      int id=std::atoi(arg[0]);
      if(strcmp(arg[1],"distance")==0){
         distance=std::atof(arg[2]);
      } else {error->all(FLERR,"Illegal transverse_pass command. Expected 'distance.'");}
      if(strcmp(arg[3],"increment")==0){
         increment=std::atof(arg[4]);
      } else {error->all(FLERR,"Illegal transverse_pass command. Expected 'increment.'");}
      transverse_passes[id]=TransversePass(distance,increment);

   } else if (strcmp(command,"am_cartesian_layer") == 0) {
      if (narg != 10) error->all(FLERR,"Illegal cartesian_layer command.");
      double x=0.0, y=0.0;
      int pass_id=-1, transverse_pass_id=-1;
      bool serpentine=false;
      int id=std::atoi(arg[0]);
      if(strcmp(arg[1],"start_position")==0){
         x=std::atof(arg[2]);
         y=std::atof(arg[3]);
      } else {error->all(FLERR,"Illegal cartesian_layer command. Expected 'start_position.'");}
      if(strcmp(arg[4],"pass_id")==0){
         pass_id=std::atoi(arg[5]);
      } else {error->all(FLERR,"Illegal cartesian_layer command. Expected 'pass_id.'");}
      if(strcmp(arg[6],"transverse_pass_id")==0){
         transverse_pass_id=std::atoi(arg[7]);
      } else {error->all(FLERR,"Illegal cartesian_layer command. Expected 'transverse_pass_id.'");}
      if(strcmp(arg[8],"serpentine")==0){
      	 //Parsing this as boolean wasn't working, temporarily changed to integer
         serpentine=std::atoi(arg[9]);
      } else {error->all(FLERR,"Illegal cartesian_layer command. Expected 'serpentine.'");}
      {
         // Create 'RectangularLayer'
         point start(x,y,0);
         Pass p=passes[pass_id];
         DIR dir=p.get_dir();
         double speed=p.get_speed();
         double pass_distance=p.get_distance();
         //Define the "overpass", which will be determined by tail_HAZ + cap_height
         double overpass = tail_HAZ + cap_HAZ;  
         TransversePass tp=transverse_passes[transverse_pass_id];
         double transverse_pass_distance=tp.get_distance();
         double transverse_pass_increment=tp.get_increment();
         rectangular_layers[id]=RectangularLayer(start,dir,speed,pass_distance,overpass,transverse_pass_distance,transverse_pass_increment,serpentine);
      }

   } else if (strcmp(command,"am_pattern") == 0) {
      int num_layers;
      vector<int> layer_ids;
      double z_start=-1.0;
      double z_increment=-1.0;
      int id=std::atoi(arg[0]);
      if(strcmp(arg[1],"num_layers")==0){
         num_layers=std::atoi(arg[2]);
      } else {error->all(FLERR,"Illegal pattern command. Expected 'num_layers.'");}
      if(strcmp(arg[3],"layer_ids")==0){
      } else {error->all(FLERR,"Illegal pattern command. Expected 'layer_ids.'");}
      int num_args=1+2+1+num_layers+2+2;
      if (narg != num_args) error->all(FLERR,"Illegal pattern command.");
      int iarg=4;
      for(int i=0;i<num_layers;i++,iarg++){
         layer_ids.push_back(std::atoi(arg[iarg]));
      }
      if(strcmp(arg[iarg],"z_start")==0){
         iarg+=1;
         z_start=std::atof(arg[iarg]);
         iarg+=1;
      } else {error->all(FLERR,"Illegal pattern command. Expected 'z_start.'");}
      if(strcmp(arg[iarg],"z_increment")==0){
         iarg+=1;
         z_increment=std::atof(arg[iarg]);
         iarg+=1;
      } else {error->all(FLERR,"Illegal pattern command. Expected 'z_increment.'");}
      pattern=Pattern(layer_ids,z_start,z_increment);

   } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAdditive::grow_app()
{
  spin = iarray[0];
  MobilityOut = darray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAdditive::init_app()
{
   delete [] sites;
   delete [] unique;
   sites = new int[1 + maxneigh];
   unique = new int[1 + maxneigh];

   dt_sweep = 1.0/maxneigh;

   int flag = 0;
   for (int i = 0; i < nlocal; i++)
    if (spin[i] < 1 || spin[i] > nspins) flag = 1;
   int flagall;
   MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
   if (flagall) error->all(FLERR,"One or more sites have invalid values");

   int next_layer_id=pattern.begin();
   active_layer=rectangular_layers[next_layer_id];

   // Compute distance function based upon initial pool position
   app_update(0.0);
}

/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */

void AppPottsAdditive::app_update(double dt)
{
   // Move pool
   if(active_layer.move(dt)){
   } else {
      // Need to update layer
      int next_layer_id=pattern.next(); 
      active_layer=rectangular_layers[next_layer_id];
   }
   // WARNING: this should always be run after checking on a move;
   // z-elevation of active layer
   double layer_z=pattern.get_layer_z_elevation();


	
	//Use the new position as input to the mobility calculation
	//Loop through all of the local sites and assign the new mobilities
	
	//Specify the shape of the melt pool and then calculate the distance at each local site.
	RASTER::pool_shape::AmEllipsoid ae(spot_width, melt_depth, melt_tail_length, cap_height, HAZ, tail_HAZ);
	
	//Go through all the local sites and calculate the distance.
   double d;
	for(int i=0;i<nlocal;i++){
			
		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};
		// Lattice point location relative to 'pool' position
		point xyz_r_p=active_layer.compute_position_relative_to_pool(XYZ,layer_z);

		//Temporary assignment of xo, xo is in the melt pool's reference frame!
		double xo[]={xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};


		if(xo[0] < 0 && xo[2] <= 0 && abs(xo[2]) <= depth_HAZ  && xo[0] > -tail_HAZ && abs(xo[1]) <= HAZ/2.0) {
	
			//If we're in the fusion zone, calculate distance
			if (abs(xo[1]) <= spot_width * 0.5 && abs(xo[0]) <= tail_HAZ) {	
				d = ae.distance(xo);
			}
			//If we're in the HAZ, calculate distance
			else if (abs(xo[1]) <= HAZ/2.0 && abs(xo[0]) <tail_HAZ) {
				d = ae.distance(xo);
			}
		}
		//If we're in front of the pool, look out to a distance cap_HAZ away
		else if (abs(xo[0]) <= cap_HAZ && xo[2] <=0 && abs(xo[2]) <= depth_HAZ && abs(xo[1]) <= HAZ/2.0) {
			d = ae.distance(xo);
		}
		else d = -10;

		//Only calculate mobilities for things inside the HAZ bounds and below the active layer
		if (d >= 0) {
			MobilityOut[i] =  compute_mobility(i, d);
		}
		//Inside the pool so set mobilty to 1 (which randomizes things)
		else if (d > -5) {

			MobilityOut[i] = 1;
		}
		//If we're outside the region of interest, make Mobility zero
		else {
			MobilityOut[i] = 0;
		}
	}
   active_layer.move(dt);
}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double AppPottsAdditive::compute_mobility(int site, double d)  {
    
	//We're going to take care of categorizing all the little details of the mobility
	//gradient in app_update, so here we'll just calculate the mobility based on distance
	MobilityOut[site] = exp(-exp_factor * d);
	
	return MobilityOut[site];
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
void AppPottsAdditive::site_event_rejection(int site, RandomPark *random) {
   int oldstate = spin[site];
   double einitial = site_energy(site);

   if (MobilityOut[site] < 0.0) {
      MobilityOut[site] = 0.0;
      return;
   }
 
   if(MobilityOut[site] >= 1.0){
      //Mobility = 0.0;
      spin[site] = (int) (nspins*random->uniform());
      return;
   }

   // events = spin flips to neighboring site different than self
   int j,m,value;
   int nevent = 0;
   int z = xyz[site][2];

   if((MobilityOut[site] > 0.0) && (MobilityOut[site] < 1.0)) {
      //(spin[i] != nspins) another criteria to exclude gg interaction
      for (j = 0; j < numneigh[site]; j++) {
         value = spin[neighbor[site][j]];
         if (value == spin[site] || value == nspins) continue;
            for (m = 0; m < nevent; m++)
               if (value == unique[m]) break;
            if (m < nevent) continue;
         unique[nevent++] = value;
      }

      if (nevent == 0) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;
         spin[site] = unique[iran];
      double efinal = site_energy(site);

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > MobilityOut[site]){
            spin[site] = oldstate;
         }
      } else if (temperature == 0.0) {
         spin[site] = oldstate;
      } else if (random->uniform() > MobilityOut[site] * exp((einitial-efinal)*t_inverse)) {
         spin[site] = oldstate;
      }

      if (spin[site] != oldstate) naccept++;

      // set mask if site could not have changed
      // if site changed, unset mask of sites with affected propensity
      // OK to change mask of ghost sites since never used

      if (Lmask) {
         if (einitial < 0.5*numneigh[site]) mask[site] = 1;
         if (spin[site] != oldstate)
         for (int j = 0; j < numneigh[site]; j++)
            mask[neighbor[site][j]] = 0;
      }
   }
}
