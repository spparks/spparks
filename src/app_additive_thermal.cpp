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
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <random>
#include "string.h"
#include "math.h"
#include "app_additive_thermal.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"
#include "comm_lattice.h"
#include "solve.h"
#include "lattice.h"
#include "timer.h"
#include "math_const.h"
#include "domain.h"
#include "output.h"

using namespace SPPARKS_NS;
using namespace MathConst;

enum{RANDOM};

/* ---------------------------------------------------------------------- */

AppAdditiveThermal::AppAdditiveThermal(SPPARKS *spk, int narg, char **arg) :
  AppPotts(spk,narg,arg)
{
  // only error check for this class, not derived classes

  double velIn;
  double x_devIn,y_devIn,z_devIn;

  if (narg != 13) error->all(FLERR,"Illegal app_style command");

  nspins = atoi(arg[1]); //Number of spins/grain IDs
  velIn = atof(arg[2]); //Velocity of travel (in meters/second)
  time_step = atof(arg[3]); //The size of each FD_step
  path_file_name = arg[4]; //Name of the input file with the proposed scan path  
  x_devIn = atof(arg[5]); //Set the standard deviation for the gaussian source. This will be in meters
  y_devIn = atof(arg[6]); //Std. dev. in y-direction
  z_devIn = atof(arg[7]); //For a volumetric Gaussian source
  short_wait_time = atof(arg[8]); //How long to pause during a layer, e.g. while laser is off and moving between rasters
  flux_prefactor = atof(arg[9]); //Total absorbed laser power in Watts
  substrate_height = atoi(arg[10]); //Substrate height in lattice sites
  nsmooth = atoi(arg[11]); //Number of Potts smoothing steps 
  dx = atof(arg[12]); //Lattice spacing (in meters)

  // read in the entire laser path

  path_file();
  done_flag =0;

  path_index = 0;
  ndouble = 3;
  allow_app_update = 1;
  ninteger = 2;
  nlocal_app = 0;
  allow_kmc = 0;
  d_residual = 0;
  vel = velIn/dx;
  x_dev = x_devIn/dx;
  y_dev = y_devIn/dx;
  z_dev = z_devIn/dx;
  recoating_time = 10;
  Kmc = 0.27695; //SPPARKS KMC scaling factor for sq_26 lattice, should not change for different materials

  //Set default values for 304L stainless steel. Can be modified in input file.
  k_solid = 30;
  k_powder = 0.3;
  boundary_temp = 300;
  density = 5706;
  T_room = 300;
  h = 25;
  eta = 1;
  Ko = 0.0000204133;
  Q = 128312;
  Tl = 1723;
  Ts = 1673;
  No = 1e15;
  Tc = 5;
  Tsig = 3;
  sizeNorm = pow(dx,3) * 2;
  sizeSig = pow(dx,3);
  latent_heat = 285000;
  specific_heat_length = 6;
  solid_front_length = 4;
  //Let's put in default values for array-based parameters
  specific_heat_temps = new double[6];
  specific_heat_vals = new double[6];
  solid_front_coeffs = new double[4];
    
  specific_heat_temps[0] = 300;
  specific_heat_temps[1] = 600;
  specific_heat_temps[2] = 1400;
  specific_heat_temps[3] = 1550;
  specific_heat_temps[4] = 1650;
  specific_heat_temps[5] = 1750;
    
  specific_heat_vals[0] = 480;
  specific_heat_vals[1] = 550;
  specific_heat_vals[2] = 675;
  specific_heat_vals[3] = 693;
  specific_heat_vals[4] = 714;
  specific_heat_vals[5] = 736;
    
  solid_front_coeffs[0] = 1.091e-5;
  solid_front_coeffs[1] = -2.034e-4;
  solid_front_coeffs[2] = 2.74e-3;
  solid_front_coeffs[3] = 1.151e-4;

  recreate_arrays();
}


/* ----------------------------------------------------------------------
   input script commands unique to this app
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"k_solid") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal solid conductivity command");
    k_solid = atof(arg[0]);
    if (k_solid <= 0) 
      error->all(FLERR,"Illegal solid conductivity");
  } 
  else if (strcmp(command,"k_powder") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal powder conductivity command");
    k_powder = atof(arg[0]);
    if (k_powder <= 0) 
      error->all(FLERR,"Illegal powder conductivity");
  }
  else if (strcmp(command,"boundary_temperature") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal boundary temperature command");
    boundary_temp = atof(arg[0]);
    if (boundary_temp <= 0) 
      error->all(FLERR,"Illegal boundary temperature");
  }
  else if (strcmp(command,"density") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal density command");
    density = atof(arg[0]);
    if (density <= 0) 
      error->all(FLERR,"Illegal density");
  }
  else if (strcmp(command,"t_room") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal room temperature command");
    T_room = atof(arg[0]);
    if (T_room <= 0) 
      error->all(FLERR,"Illegal room temperature");
  }
  else if (strcmp(command,"emissivity") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal emissivity command");
    eta = atof(arg[0]);
    if (eta <= 0) 
      error->all(FLERR,"Illegal emissivity");
  }
  else if (strcmp(command,"convection_coefficient") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal convection coefficient command");
    h = atof(arg[0]);
    if (h <= 0) 
      error->all(FLERR,"Illegal convection coefficient");
  }
  else if (strcmp(command,"Ko") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal Arrhenius pre-factor command");
    Ko = atof(arg[0]);
    if (eta <= 0) 
      error->all(FLERR,"Illegal Arrhenius pre-factor");
  }
  else if (strcmp(command,"Q") == 0) {
    if (narg != 1) 
      error->all(FLERR,"Illegal Arrhenius exponential-factor command");
    Q = atof(arg[0]);
    if (Q <= 0) 
      error->all(FLERR,"Illegal Q");
  }
  else if (strcmp(command,"liquidus") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal liquidus command");
    Tl = atof(arg[0]);
    if (Tl <= 0) 
      error->all(FLERR,"Illegal liquidus temperature");
  }
  else if (strcmp(command,"solidus") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal solidus command");
    Ts = atof(arg[0]);
    if (Ts <= 0) 
      error->all(FLERR,"Illegal solidus temperature");
  }
  else if (strcmp(command,"nucleation_density") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal nucleation density command");
    No = atof(arg[0]);
    if (No < 0) 
      error->all(FLERR,"Illegal nucleation density");
  }
  else if (strcmp(command,"critical_undercooling") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal critical_undercooling command");
    Tc = atof(arg[0]);
    if (Tc < 0) 
      error->all(FLERR,"Illegal critical undercooling");
  }
  else if (strcmp(command,"undercooling_deviation") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal undercooling_deviation command");
    Tsig = atof(arg[0]);
    if (Tsig < 0) 
      error->all(FLERR,"Illegal undercooling standard deviation");
  }
  else if (strcmp(command,"mean_nuclei_volume") == 0) {
    if (narg !=1) error->all(FLERR,"Illegal mean_nuclei_volume command");
    sizeNorm = atof(arg[0]);
    if (sizeNorm < 0) 
      error->all(FLERR,"Illegal mean nuclei volume");
  }
  else if (strcmp(command,"nuclei_volume_deviation") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal nuclei_volume_deviation command");
    sizeSig = atof(arg[0]);
    if (sizeSig < 0) 
      error->all(FLERR,"Illegal nuclei volume standard deviation");
  

  } else if (strcmp(command,"specific_heat") == 0) {
    delete [] specific_heat_temps;
    delete [] specific_heat_vals;
    if (narg < 1) error->all(FLERR,"Illegal specific_heat command");
    int nHeat;
    nHeat = atoi(arg[0]);
    if (nHeat <= 0) error->all(FLERR,"Illegal specific heat specification");
    specific_heat_temps = new double [nHeat];
    specific_heat_vals = new double [nHeat];
    int j = 0;
    for(int i = 1; i < nHeat * 2 + 1; i = i + 2) {
      specific_heat_temps[j] = atof(arg[i]);
      specific_heat_vals[j] = atof(arg[i + 1]);
      j++;
    }

  } else if (strcmp(command,"latent_heat") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal latent_heat command");
    latent_heat = atof(arg[0]);
    if (latent_heat < 0) error->all(FLERR,"Illegal latent heat value");
    
    // need to modify to get solid_front_length

  } else if (strcmp(command,"solid_front_vel") == 0) {
    delete [] solid_front_coeffs;
    if (narg < 1) error->all(FLERR,"Illegal solid_front_vel command");
    solid_front_length = atoi(arg[0]);
    if (solid_front_length <= 0)
      error->all(FLERR,"Illegal solidification front velocity specification");
    solid_front_coeffs = new double [solid_front_length];
    int j = 0;
    for (int i = 1; i < solid_front_length + 1; i++) {
      solid_front_coeffs[j] = atof(arg[i]);
      j++;
    }

  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::grow_app()
{
  spin = iarray[0];
  activeFlag = iarray[1];
  MobilityOut = darray[0];
  T = darray[1];
  SolidD = darray[2];
  
  if (nlocal_app < nlocal) nlocal_app = nlocal;          
}


/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::init_app()
{
  delete [] sites;
  delete [] unique;
  double sqrt2 = 1.4142135624;
  double sqrt3 = 1.7320508076;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];
  nucleationFlags = new int[nspins];
  nucleationTemps = new double[nspins];
  nucleationSizes = new double[nspins];
  
  dt_sweep = 1.0/maxneigh;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (spin[i] < 1 || spin[i] > nspins) {
      flag = 1;
    }
    //If we're above the substrate height, randomize the spins
    if(xyz[i][2] > substrate_height) {
      spin[i] = (int) (nspins*ranapp->uniform());
      activeFlag[i] = 0;
    }
    //If we're less than the value, set activeFlag to the "solid" condition
    else {
      activeFlag[i] = 3;
    }
  }
  comm->all();
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
  
  //Initialize the nucleationFlags vector
  if (domain->me==0) {
    nucleation_spins(ranapp);    
  }
  
  MPI_Bcast(nucleationFlags,nspins, MPI_INT,0,world);
  
  //Initialize the nucleationTemps and nucleationSizes vectors
  if (domain->me==0) {
    nucleation_init();    
  }
  
  MPI_Bcast(nucleationTemps,nspins, MPI_DOUBLE,0,world);
  MPI_Bcast(nucleationSizes,nspins, MPI_DOUBLE,0,world);
  
  //Initialize the neighDist array need to fill with good values
  neighDist = new double[26];
  neighDist[0] = sqrt3 * dx;
  neighDist[1] = sqrt2 * dx;
  neighDist[2] = sqrt3 * dx;
  neighDist[3] = sqrt2 * dx;
  neighDist[4] = dx;
  neighDist[5] = sqrt2 * dx;
  neighDist[6] = sqrt3 * dx;
  neighDist[7] = sqrt2 * dx;
  neighDist[8] = sqrt3 * dx;
  neighDist[9] =  sqrt2 * dx;
  neighDist[10] = dx;
  neighDist[11] =  sqrt2 * dx;
  neighDist[12] = sqrt3 * dx;
  neighDist[13] = sqrt3 * dx;
  neighDist[14] =  sqrt2 * dx;
  neighDist[15] = dx;
  neighDist[16] =  sqrt2 * dx;
  neighDist[17] = sqrt3 * dx;
  neighDist[18] = sqrt2 * dx;
  neighDist[19] = sqrt3 * dx;
  neighDist[20] =  sqrt2 * dx;
  neighDist[21] = dx;
  neighDist[22] =  sqrt2 * dx;
  neighDist[23] = sqrt3 * dx;
  neighDist[24] = sqrt2 * dx;
  neighDist[25] = sqrt3 * dx;
	
  //Check that our timestep is short enough to capture solidification behavior

  double max_front_vel = 0;
  int power = solid_front_length -1;
  for (int k = 0; k < solid_front_length; k++) {
    max_front_vel = max_front_vel + solid_front_coeffs[k] * pow(Tl - Ts, power);
    power--;
  }
  if (max_front_vel * time_step > dx) {
    double max_time_step = dx/max_front_vel;
    fprintf(screen,"time_step is too large to capture moving solidification front. Max allowable timestep is %f\n",max_time_step);
    error->all(FLERR,"Illegal time_step");
  }
  
  this->app_update(0.0);
}

/* ----------------------------------------------------------------------
   Read in the coordinate file before going into MPI
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::path_file()
{
  int x_max = 0;
  int y_max = 0;
  int z_max = 0;
   
  std::string line;
  std::ifstream myfile(path_file_name);
   
  //Check if the file exists and is readable.
  if (myfile.good()) {
    if(domain->me == 0) {
      fprintf(screen,"Path input file found \n");
    }
  }
  else {
    error->all(FLERR,"Path input file not found");;
  }

  // new lines will be skipped unless we stop it from happening:    
  myfile.unsetf(std::ios_base::skipws);

  int aNumOfLines = 0;

  std::string aLineStr;
  while (getline(myfile, aLineStr))
    {
      if (!aLineStr.empty())
        line_count++;
    }
	
  //Lets just use a 1D vector with fancy indexing to do our scan_array
  //Each point will be accessed by scan_array[y * sizeX + x], where x is 0-3 and y is 0 - line_count - 1
  //We should also make new arrays for each xyz and value, cause I want to easily access the variables.
  //We might not need d_scan_array, but its also not bad to have...
  scan_array = new double[line_count * 5];
  x_scan_array = new double[line_count];
  y_scan_array = new double[line_count];
  z_scan_array = new double[line_count];
  d_scan_array = new double[line_count];
  p_scan_array = new double[line_count];
	
  //Read in the input file. It should always have 5 columns X,Y,Z distance, and pause_flag
  std::ifstream in_file(path_file_name);
  
  for (int row = 0; row < line_count; row++) {
  
    getline(in_file, line);
    if(!in_file.good() )
      break;
  	
    std::stringstream iss(line);
  	
    //if(domain->me == 0) {
    //	std::cout << line << "\n";
    //}

  	
    for (int col = 0; col < 5; col++) {
  		
      std::string val;
      getline(iss, val, ',');

  		
      std::stringstream convertor(val);

      convertor >> scan_array[row * 5 + col];  		
    }
  }
  
  for (int i = 0; i < line_count; i++) {
  
    x_scan_array[i] = scan_array[i * 5 + 0];
    y_scan_array[i] = scan_array[i * 5 + 1];
    z_scan_array[i] = scan_array[i * 5 + 2];
    d_scan_array[i] = scan_array[i * 5 + 3];
    p_scan_array[i] = scan_array[i * 5 + 4];

  }
  	
  delete [] scan_array;
  in_file.close();
  
}

/* ----------------------------------------------------------------------
   perform finite difference on a single site
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::site_event_finitedifference(int i)
{
    
  //set values used for convenience
  double chm_temp = 0.0;
	
  int good_neigh [ ] = {4, 10, 12, 15, 21, 13};
  int oppo_neigh [ ] = {21, 15, 13, 10, 4, 12};
	
  double local_temp = 0.0;
  int num_inactive = 0;
  double Cp = specificHeatCalculator(T[i]);
  double sigma = 5.670367e-8;
  double sigma_adj = eta * sigma;
  double flux_loc;
	
  //Figure out the flux value for the current heat source.
  flux_loc = flux_finder(i);
	
  //Check if we're on a global boundary and fix the temperature
  //Kyle uses convection boundary conditions at the bottom of his substrate. Maybe something we should add or change to.
  if(xyz[i][0] == domain->boxxlo || xyz[i][0] > domain->boxxhi - 2 || xyz[i][1] == domain->boxylo || xyz[i][1] > domain->boxyhi -2 || xyz[i][2] == domain->boxzlo) {
    T[i] = boundary_temp;
    return;
  }
	
  //First thing, loop through the neighbors of a site and see if any of them are inactive
  //We need to deal with when a site has more than one inactive value, which isn't that uncommon...
  //Lets average the contributions from all the orientations
  for (int j = 0; j < 6; j++) {
		
    int s = neighbor[i][good_neigh[j]];
		
    //Loop through neighbors, do usual FD if active, and do BCs if not
    //Treat solid & liquid cases in the same way
    //Site is Solid (or liquid)
    if( activeFlag[i] == 3 || activeFlag[i] == 2) {
      //Solid-solid
      if(activeFlag[s] == 3 || activeFlag[s] == 2) {
        chm_temp += k_solid/(Cp * density)  * (T[s] - T[i])/(dx*dx);
      }
      //Solid-powder, average conductivities
      else if (activeFlag[s] == 1) {
        chm_temp += (k_solid + k_powder)/(2*Cp*density) * (T[s] - T[i])/(dx*dx);
      }
      //Do convection-radiation
      else if (activeFlag[s] == 0) {
        //Implement Convection BCs
        num_inactive++;
        int oppo = neighbor[i][oppo_neigh[j]];
        double Tb = (3 * T[i] - T[oppo])/2;
        chm_temp += (h * (Tb - T_room) + eta * sigma * (Tb * Tb * Tb * Tb - T_room * T_room * T_room * T_room))/density/Cp/dx;
      }
    }
    //Site is powder
    else if(activeFlag[i] == 1) {
      //powder-powder
      if(activeFlag[s] == 1) {
        chm_temp += k_powder/(Cp * density)  * (T[s] - T[i])/(dx*dx);
      }
      //Powder-solid, average conductivities
      else if (activeFlag[s] == 3 || activeFlag[s] == 2) {
        chm_temp += (k_solid + k_powder)/(2*Cp * density) * (T[s] - T[i])/(dx*dx);
      }
      //Do convection-radiation
      else if (activeFlag[s] == 0) {
        //Implement Convection BCs
        num_inactive++;
        int oppo = neighbor[i][oppo_neigh[j]];
                
        double Tb = (3 * T[i] - T[oppo])/2;
        chm_temp += (h * (Tb - T_room) + eta * sigma * (Tb * Tb * Tb * Tb - T_room * T_room * T_room * T_room))/density/Cp/dx;
      }
    }
  }
	
  //Add the laser power/volume
  chm_temp += flux_loc/(Cp * density);
			
  T[i] += time_step * chm_temp;
  return;
	
}

//Represent specific heat as a piecewise function, so we can account for the rapid increase
//upon melting (corresponding to release of latent heat) and its larger value in the liquid phase.
//This function will take a temperature value as input and return the corresponding specific heat value

double AppAdditiveThermal::specificHeatCalculator(double temper) {

  double slope = 0;
  double Cp = 0;
    
  //Check if we are above or below the specified temperature range. If not, linearly interpolate the proper value.
  if(temper >= specific_heat_temps[specific_heat_length -1]) {
    Cp = specific_heat_vals[specific_heat_length -1];
  }
  else if(temper <= specific_heat_temps[0]) {
    Cp = specific_heat_vals[0];
  }
  else {
    for(int i = 0; i < specific_heat_length; i++) {
      if(i == 0 && temper <= specific_heat_temps[0]) {
        Cp = specific_heat_vals[0];
        break;
      }
      else if (temper > specific_heat_temps[specific_heat_length -1]) {
        Cp = specific_heat_vals[i];
        break;
      }
      else if (temper > specific_heat_temps[i -1] && temper <= specific_heat_temps[i]) {
        slope = (specific_heat_vals[i] - specific_heat_vals[i -1])/(specific_heat_temps[i] - specific_heat_temps[i-1]);
        Cp = slope * (temper - specific_heat_temps[i-1]) + specific_heat_vals[i -1];
        break;
      }
    }
    //Add latent heat if in mushy zone
    if (temper < Tl && temper > Ts) {
      Cp += latent_heat/(Tl - Ts);
    }
  }
    
  return Cp;
}


/* ----------------------------------------------------------------------
   At each timestep, we will go through and calculate the temperature field from
   the total distance traveled. We need to account for negative numbers, and 
   beginning/end conditions. At any point if the value is negative, we do everything
   but then just throw the mobility values away.
   We always measure local distance from the point with the smallest index.
   We're currently starting from the beginning of the array everytime we search for our location, which is dumb.
   What would be better is to either use a linked list/stack where we can remove
   the index after using it, or keep track of the index of the current start of our
   interval and go from there.
   We should be going a constant dt in time between each step.
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::position_finder_in() {
  double local_dist;
  double x_vec;
  double y_vec;
  int longer_index = 0;
	
  //We should have a special case for when we restart a simulation during a "dwell" period
  if(path_index == 0 && p_scan_array[0] != 0) {
		
    //I think we should just restart running the entire wait time for now.
    //Its not too simple with our current set up to adjust the amount of time left.
    if (p_scan_array[0] < 1.5) {
      wait_time = recoating_time;
      if (domain->me == 0) fprintf(screen, "Recoating\n");
    }
    else {
      wait_time = short_wait_time;
      if (domain->me == 0) fprintf(screen, "Moving laser\n");
    }
    wait_time = wait_time - time_step;
    x_meltspot = domain->boxxhi * 1000;
    y_meltspot = 0;
    z_meltspot = z_scan_array[path_index];
    return;
  }


  //First check if in a pause block and have time to  wait
  if(p_scan_array[path_index] != 0 && wait_time > 0) {
		
    wait_time = wait_time - time_step;
    x_meltspot = domain->boxxhi * 1000;
    y_meltspot = 0;
    z_meltspot = z_scan_array[path_index];
    return;
  }
  //If we've been waiting but shouldn't wait anymore, start at next index
  else if (p_scan_array[path_index] != 0) {
    path_index = path_index + 1;
    wait_time = 0;
    d_residual = 0;
		
    x_meltspot = x_scan_array[path_index];
    y_meltspot = y_scan_array[path_index];
    z_meltspot = z_scan_array[path_index];
    return;
  }
  //We're not starting in a wait spot, need to loop through until we find the active interval
  else {
    //Just advance by vel everytime this is called.
    double long_local_dist = d_residual + vel * time_step + d_scan_array[path_index];
		
    for(int i = 1; i < line_count - path_index; i++) {
			
      //If we've crossed into a pause area, stop and interval beginning and set flags
      if(p_scan_array[path_index + i - 1] != 0) {
        path_index = path_index + i -1;
        if (p_scan_array[path_index] < 1.5) {
          wait_time = recoating_time;
          if (domain->me == 0) fprintf(screen, "Recoating\n");
        }
        else {
          wait_time = short_wait_time;
          if (domain->me == 0) fprintf(screen, "Moving laser\n");
        }
        d_residual = 0;
				
        x_meltspot = x_scan_array[path_index];
        y_meltspot = y_scan_array[path_index];
        z_meltspot = z_scan_array[path_index];	
        return;
      }
			
      //Otherwise, if we've found an index with a longer distance than our current one, do the previous thing
      else if(d_scan_array[path_index + i] > long_local_dist)  {
        longer_index = path_index + i - 1;
        local_dist = long_local_dist - d_scan_array[path_index + i - 1];
        path_index = path_index + i - 1;
        d_residual = local_dist;
        wait_time = 0;
        break;
      }
      else if(line_count == path_index + 2){
        done_flag = 1;
        x_meltspot = x_scan_array[path_index];
        y_meltspot = y_scan_array[path_index];
        z_meltspot = z_scan_array[path_index];
        return;
      }
    }
  }

  //We'll find the vector between the two end points,
  //use this to find the point's distance from the plane, and then test if its within the ellipsoid.
  x_vec = abs(x_scan_array[longer_index + 1]) - abs(x_scan_array[longer_index]);
  y_vec = y_scan_array[longer_index + 1] - y_scan_array[longer_index];
 	
  //We can calculate the variables on the spine quite easily by taking the normalized vector, multiplying it by the small length, and adding it to the endpoint.
  x_meltspot = x_vec/(d_scan_array[longer_index + 1] - d_scan_array[longer_index]) * local_dist + abs(x_scan_array[longer_index]);
  y_meltspot = y_vec/(d_scan_array[longer_index + 1] - d_scan_array[longer_index]) * local_dist + abs(y_scan_array[longer_index]);
  z_meltspot = z_scan_array[longer_index];
  return;
}

/* ----------------------------------------------------------------------
   Compute the mobility at the specified lattice site. Returns a double
   between 0 and 1 representing the mobility. This is used as normal to compute
   the mobility during solid-state grain growth.
   The MobilityOut array that the values are assigned to is also re-used to control
   smoothing steps immediately after solidification by using negative values
   ------------------------------------------------------------------------- */

double AppAdditiveThermal::compute_mobility(int i, RandomPark *random)  {
    
  //Don't assign the global array in this function
  double mob;
	 
  //If we're a solid phase, set mobility for grain growth
  if (activeFlag[i] == 3) {
    mob = exp(-Q/(R*T[i]));
  }
  else {
    error->all(FLERR,"Error: compute_mobility called on non-eligible site. Finite difference solver is likely unstable. Try decreasing timestep.");
  }
  return mob; 
}   

/* ----------------------------------------------------------------------
   Simple function to calculate flux at a given lattice site.
   This should represent the flux going into an element the size of of a site
   ------------------------------------------------------------------------- */

double AppAdditiveThermal::flux_finder(int site) {
	
  double flux_loc;
  double preFactorTrue = 2 * flux_prefactor/((x_dev * dx * MY_2PIS) * (y_dev * dx * MY_2PIS) * (z_dev * dx * MY_2PIS));
	
  //Normalized gaussian. Total flux will stay constant no matter the values of x,y,z_dev
  flux_loc =  preFactorTrue * exp(-0.5 * ((pow((xyz[site][0] - x_meltspot), 2)/(x_dev*x_dev) + pow((xyz[site][1] - y_meltspot), 2)/(y_dev*y_dev) + pow((xyz[site][2] - z_meltspot),2)/(z_dev * z_dev))));
	
  //If we're above the current layer or have reached the end of the path file, set the flux to zero.
  if (xyz[site][2] > floor(z_meltspot) || done_flag == 1) {
    flux_loc = 0;
  }

  return flux_loc;
}

/* ----------------------------------------------------------------------
   iterate through the temperature solver and update phases as needed
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::app_update(double dt)
{
  double localtime=0.0;
  double tempMax = 0.0;
  double tempMaxAll = 0.0;

  //communicate all sites to make sure it's up-to-date when it starts
  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);

  //Calculate the current position
  if(domain->me == 0) {
    position_finder_in();
    //         fprintf(screen,"path_index %d, line_count %d\n", path_index, line_count);
  }

  timer->stamp(TIME_APP);
  MPI_Bcast(&done_flag,1,MPI_INT,0,world);
  MPI_Bcast(&path_index,1,MPI_INT,0,world);
  MPI_Bcast(&x_meltspot,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&y_meltspot,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&z_meltspot,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&wait_time,1,MPI_DOUBLE,0,world);
  timer->stamp(TIME_COMM);

  //If we're recoating, run until things get below Ts and then reset all values to room temp
  if(wait_time > 1) {
    //See if we're below melting throughout the domain
    tempMax = compute_tempMax();
    MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
    //If below, set all values to room temp and skip to the end of the wait
    if(tempMaxAll < Ts) {
      if(domain->me == 0) {
        fprintf(screen,"Continuing after waiting %f\n",recoating_time - wait_time);
      }
      wait_time = -1;
      for (int i=0; i<nlocal; i++) {
        T[i] = boundary_temp;
      }
      comm->all();
      return;
    }
  }

  for (int i=0; i<nlocal; i++) {
    //Keep looping if we're above the meltspot
    if(xyz[i][2] > floor(z_meltspot)) continue;
        
    //If below melt spot, run finite difference
    //This could be slightly screwy for the first step after a new layer
    else if (activeFlag[i] == 0) {
      activeFlag[i] = 1;
    }
    site_event_finitedifference(i);
        
    //Let's also update the active flag after each FD loop
    //This is also the place to handle nucleation and solidification front impingement
    if(activeFlag[i] == 3) {
      if(T[i] > Tl) {
        activeFlag[i] = 2;
        //Also randomize spin and reset cumulative variables
        spin[i] = (int) (nspins * ranapp->uniform());
        SolidD[i] = 0;
        MobilityOut[i] = 0;
      }
    }
    //If we're molten, call the mushy_phase function to figure out any phase change
    else if (activeFlag[i] == 2 && T[i] <= Tl) {
      mushy_phase(i, ranapp);
    }
    //Go from powder to molten
    else if (activeFlag[i] == 1 && T[i] > Tl) {
      activeFlag[i] = 2;
      //Also randomize spin and reset cumulative variables
      spin[i] = (int) (nspins * ranapp->uniform());
      SolidD[i] = 0.0;
      MobilityOut[i] = 0;
    }
  }
  timer->stamp(TIME_APP);

  //re-sync all the data
  comm->all();
  timer->stamp(TIME_COMM);

	
  //Check if this is the last iteration, if so, output an updated version of the path file.
  if(time + dt >= stoptime && domain->me == 0) {	
    path_file_update();
  }
	
}
/* ----------------------------------------------------------------------
   We'll often need to restart these simulations. We calculate our spot location by tracking
   an index. Therefore, it'll be easiest to rewrite the path file so that it restarts at our
   final location.
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::path_file_update()
{
  //Output a restart_path file.
  std::string restart_file_str("restart_path.txt");
  FILE* fout = fopen(restart_file_str.c_str(), "w");
	
  double start_distance = d_scan_array[path_index] + d_residual;
	
  //The first index value should be our current location, not array value
  //We might be off by one for the first index
  if(p_scan_array[path_index] == 0) {
    x_scan_array[path_index] = x_meltspot;
    y_scan_array[path_index] = y_meltspot;
    z_scan_array[path_index] = z_meltspot;
  }	
  d_scan_array[path_index] = 0;
	
	
  for( int i = path_index; i < line_count; i++) {
    //Update the scan distance
    if(i > path_index) {
      d_scan_array[i] = d_scan_array[i] - start_distance;
    }
				
    fprintf(fout, "%f,\t%f,\t%f,\t%f,\t%f\n", x_scan_array[i],y_scan_array[i],z_scan_array[i],d_scan_array[i],p_scan_array[i]);
  }
  fclose(fout);
}


/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy(i);
  double Mobloc = 0;
    
  //Assign the local mobility
  Mobloc = MobilityOut[i];

    
  int j,m,value;
  int nevent = 0;
  
  if((Mobloc < 0.0) || (Mobloc > 1.0001)) {
    MobilityOut[i] = 0;
    return;
  }


  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    //Exclude gas, powder or molten sites from the Potts neighbor tally
    if (value == spin[i] || value == nspins || (activeFlag[neighbor[i][j]] != 3 && activeFlag[neighbor[i][j]] != 1)) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform());
  if (iran >= nevent) iran = nevent-1;
  spin[i] = unique[iran];
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
    if (random->uniform() > Mobloc){
      spin[i] = oldstate;
    }
  }
  else if (temperature == 0.0) {
    spin[i] = oldstate;
  } 
  else if (random->uniform() > Mobloc * exp((einitial-efinal)*t_inverse)) {
    spin[i] = oldstate;
  }

  if (spin[i] != oldstate) naccept++;

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (spin[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
  	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   Perform evolution for sites in the mushy zone. There are several things that need to happen.
   1. Determine if the site is solid or liquid (from activeFlag)
   2. Determine if the site should nucleate a new grain (from nucleationFlags)
   3. If our current site is liquid & can't nucleate, have it try to switch to a solid neighbor
   (with 4 calculated from the undercooling in someway, not temperature)
   4. If our current site is liquid & can nucleate, check if local undercooling is equal to its 
   critical temp. If so, change the activeFlag value to solid. If not, see if there are any solid sites
   that should capture it.
   5. If our current site is solid, see if it should flip to a neighboring solid value (with
   mobility calculated from undercooling.)
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::mushy_phase(int i, RandomPark *random){
  int nevent = 0;
  int m,value;
  double Tcool = Tl - T[i];
  //For default settings, SolidD[i] =+ (1.091e-5 * pow(Tcool, 3) - 2.034e-4 * pow(Tcool, 2) + 2.74e-3*Tcool + 1.151e-4) * time_step;
       
  //Our site should always be molten and below Tl
  //Check if it's eligible to nucleate
  if(nucleationFlags[spin[i]]) {
    //Can and will nucleate
    if(Tcool >= nucleationTemps[spin[i]]){
      activeFlag[i] = 3;
      //Don't let nucleated site disappear during smoothing
      SolidD[i] = -nsmooth-2;
      return;
    }
    //Can nucleate, but won't yet. Allow to if the solidification front gets captured.
    else {
      //Add the distance of the front travel. This is for 304L. Need to multiply by timestep to get distance
      //Try doing this with an arbitrary array
      int power = solid_front_length -1;
      for(int k = 0; k < solid_front_length; k++) {
        SolidD[i] = SolidD[i] + solid_front_coeffs[k] * pow(Tcool, power) * time_step;
        power--;
      }
      //Go through neighbor list and add them to possible switches
      for (int j = 0; j < numneigh[i]; j++) {
        if(neighDist[2] <= SolidD[i] && activeFlag[neighbor[i][j]] == 3) {
          value = spin[neighbor[i][j]];
          //Exclude gas or molten sites from the Potts neighbor tally
          if(neighDist[j] > SolidD[i]) continue;
          unique[nevent++] = value;										
        }
      }
      //If no neighbor is eligible, return before changing anything. Will try next sweep.
      if (nevent == 0) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;
      spin[i] = unique[iran];
      activeFlag[i] = 3;
      SolidD[i] = -1;
      return;
    }
  }
  else {
    //Add the distance of the front travel. This is for 304L I think...
    int power = solid_front_length -1;
    for(int k = 0; k < solid_front_length; k++) {
      SolidD[i] = SolidD[i] + solid_front_coeffs[k] * pow(Tcool, power) * time_step;
      power--;
    }
       
    //Go through neighbor list and add them to possible switches
    for (int j = 0; j < numneigh[i]; j++) {
      if(neighDist[2] <= SolidD[i] && activeFlag[neighbor[i][j]] == 3) {
        value = spin[neighbor[i][j]];
        if(neighDist[j] > SolidD[i]) continue;
        unique[nevent++] = value;										
      }
    }
        
    //If no neighbor is eligible, return before changing anything. Will try next sweep.
    if (nevent == 0) return;
    int iran = (int) (nevent*random->uniform());
    if (iran >= nevent) iran = nevent-1;
    spin[i] = unique[iran];
    activeFlag[i] = 3;
    SolidD[i] = -1;
    return;
  }
}

/* ----------------------------------------------------------------------
   Only nucleating one site at a time introduce lattice size dependency. Here we will
   use a user-defined nucleation particle size and flip neighboring sites until that size is met
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::nucleation_particle_flipper(int i, int partRad, RandomPark *random) {
    
  //If one site is big enough to satisfy, skip evertyhing
  if(partRad <= 0) return;
    
  int nSites = partRad;
  int nSitesIn = nSites;
  int nearest_neigh [ ] = {4, 10, 12, 15, 21, 13};
  int second_nearest [ ] = {9,16,3,22,14,11,20,5,1,24,7,18};
  int third_nearest [ ] = {0,25,23,2,6,19,17,8};
  int nneigh = 0;
  int possible_neigh[26];
    
  //Its hard to go through shells iteravely if the number of neighbors isn't the full 26.
  //Check if this is the case and just loop through the neihgbor list if so
  if(numneigh[i] != 26) {        
    for(int j = 0; j < numneigh[i]; j++) {
      if(activeFlag[neighbor[i][j]] == 2) {
        spin[neighbor[i][j]] = spin[i];
        activeFlag[neighbor[i][j]] = 3;
        SolidD[neighbor[i][j]] = -nsmooth -3;
        nSites--;
      }
      if(nSites <= 0) {
        return;
      }
    }
  }
  else {
    //Go through neighbors and nucleate liquid ones
    //Do 1st ones first (this isn't random yet)
    for(int j = 0; j < 6; j++) {
      if(activeFlag[neighbor[i][nearest_neigh[j]]] == 2) {
        spin[neighbor[i][nearest_neigh[j]]] = spin[i];
        activeFlag[neighbor[i][nearest_neigh[j]]] = 3;
        SolidD[neighbor[i][nearest_neigh[j]]] = -nsmooth -3;
        nSites--;
      }
      if(nSites <= 0) {
        return;
      }
    }
    //Do 2nd shell
    for(int j = 0; j < 12; j++) {
      if(activeFlag[neighbor[i][second_nearest[j]]] == 2) {
        spin[neighbor[i][second_nearest[j]]] = spin[i];
        activeFlag[neighbor[i][second_nearest[j]]] = 3;
        SolidD[neighbor[i][second_nearest[j]]] = -nsmooth -3;
        nSites--;
      }
      if(nSites <= 0) {
        return;
      }      
    }
    //Do 3rd shell    
    for(int j = 0; j < 8; j++) {
      if(activeFlag[neighbor[i][third_nearest[j]]] ==2) {
        spin[neighbor[i][third_nearest[j]]] = spin[i];
        activeFlag[neighbor[i][third_nearest[j]]] = 3;
        SolidD[neighbor[i][third_nearest[j]]] = -nsmooth -3;
        nSites--;
      }
      if(nSites <= 0) {
        return;
      }
    }
  }
  //If we didn't fill any sites this time, just quit
  if (nSites == nSitesIn) {
    return;
  }
    
  //If we still haven't satisfied the particle size, pick a neighbor at random and solidify from there.
  //Build a list of same-particle neighbors and pick one randomly
  for(int j =0; j < numneigh[i]; j++) {
    if(spin[neighbor[i][j]] == spin[i] && activeFlag[neighbor[i][j]] == 3) {
      possible_neigh[nneigh] = j;
      nneigh++;
    }
  }
  //If no possible nieghbors, quit
  if(nneigh == 0) {
    return;
  }
  //Otherwise, randomly pick a possilbe neighbor
  int neighran =  round(((nneigh -1) * random->uniform()));
  nucleation_particle_flipper(neighbor[i][possible_neigh[neighran]],nSites, random);
  return;    
}


/* ----------------------------------------------------------------------
   Nucleation site initializer. Find the volume of a voxel from dx^3 and Multiply by No.
   This will be the average number of nucleation sites in the voxel. This is also the
   fraction of spins that we want to be able to nucleate new grains. If the value is greater
   than 1 (which we should avoid), allow all spins to nucleate. If not, call a random number
   between zero and one. If the number is less than the fraction, make true. If not, make false.
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::nucleation_spins(RandomPark *random) {
  double nucleationFraction = dx * dx * dx * No;
    
  //Make all spins nucleation sites. Should avoid this.
  if(nucleationFraction >= 1.0) {
    fprintf(screen,"Nucleation fraction (%f) is greater than 1. Decrease No or increase mesh resolution.\n", nucleationFraction);
    for (int i = 0; i < nspins; i++) {
      nucleationFlags[i] = 1;
    }
  }
  //Do a random number test and allow the spin to nucleate if less than
  else {
    fprintf(screen,"Nucleation Fraction is %f \n", nucleationFraction);
    for (int i = 0; i < nspins; i++) {
      if(random->uniform() <= nucleationFraction) {
        nucleationFlags[i] = 1;
      }
      else {
        nucleationFlags[i] = 0;
      }
    }
  }   
}

/* ----------------------------------------------------------------------
   The first version of this just initialized critical nucleation temperatures.
   This version will also initialize nucleii size (starting with a normal dist)
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::nucleation_init() {
  std::normal_distribution<> dist_T{Tc,Tsig};
  std::normal_distribution<> dist_S{sizeNorm,sizeSig};
  std::random_device rd{};
  std::mt19937 gen{rd()};
    
  //Randomly assign a temperature to every spin
  for(int i = 0; i < nspins; i++) {
    nucleationTemps[i] = dist_T(gen);
    nucleationSizes[i] = dist_S(gen);
  }
}

/* ----------------------------------------------------------------------
   Evolve simulation time including rejection KMC solver. This function controls
   the relative occurence of solidification/thermal timesteps and solid-state kMC steps
   by comparing the dtMC and time_step variables. time_step is a constant variable,
   while dtMC depends on lattice size and maximum solid-phase temperature.
   ------------------------------------------------------------------------- */

void AppAdditiveThermal::iterate_rejection(double stoptime)
{
  int i,icolor,nselect,nrange,jset, nMC;
  int *site2i;
  double tempMax;
  double tempMaxAll;
  double FDElapsed;
  double nucVolume;
  double mobMax = 0;
  
  // set loop is over:
  // sectors if there are sectors and no colors
  // colors if there are colors and no sectors
  // first nsector sets if there are both sectors and colors

  int nset_loop = nset;
  if (bothflag) nset_loop = nsector;

  int done = 0;
  //This while loop is for the entire simulation run time! (stoptime = "run stoptime")
  while (!done) {
    //Find the highest temperature in the local array and compare with others
    tempMax = compute_tempMax();
    timer->stamp(TIME_APP);
    MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
    timer->stamp(TIME_COMM);
    //Using the global maximum temperature, compute our smallest timestep
    dtMC = compute_timeMin(tempMaxAll);
    //Also using this value, compute the maximum solid-state Mobility
    mobMax = exp(-Q/(R*tempMaxAll));
  	  	
    //Figure out how many loops we should do for each step
    //If Monte Carlo time is bigger than therm, do app_update loops first and do one MC step
    if(dtMC > time_step) {
  		
      FDElapsed = 0;
      nMC = 1;
      //Zero out the mobility
      for(int i = 0; i < nlocal; i++) {
        MobilityOut[i] = 0;
      }
      mobMax = 0;
      //         if(domain->me == 0) {
      //             fprintf(screen,"Potts time_step %.5e is bigger than solidification time_step %.5e.\n",dtMC, time_step);
      //         }  		
      while(dtMC >= FDElapsed) {
        //Update temperatures and phases
        app_update(time_step);
        //Find the new dtMC value
        tempMax = compute_tempMax();
        timer->stamp(TIME_SOLVE);
        MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
        timer->stamp(TIME_COMM);
        //Using the global maximum temperature, compute our smallest timestep
        //Only update dtMC if it is smaller (higher temperature)
        //Basing MC calculation off of the highest temp observed in time_step
        if (dtMC > compute_timeMin(tempMaxAll) || mobMax < 1e-8) {
          dtMC = compute_timeMin(tempMaxAll);
          //If new mobMax is larger, use it. This is tested by the enclosing if statement
          mobMax = exp(-Q/(R*tempMaxAll));
        }

        //Add to running mobility values at each site. Multiply mobility by timestep_size
        //which makes the integral of a constant function
        for(int i = 0; i < nlocal; i++) {
          if(activeFlag[i] < 3) continue;
          MobilityOut[i] += time_step * compute_mobility(i,ranapp);
          //Also check if we're just solidified and should be "relaxed"
          if(SolidD[i] < 0 && SolidD[i] > -nsmooth -1)    {
            MobilityOut[i] = 1;
            site_event_rejection(i, ranapp);
            SolidD[i]--;
          }
          //Check if we just nucleated and need to grow larger
          else if(SolidD[i] == -nsmooth -2) {
            nucVolume = nucleationSizes[spin[i]];
            nucleation_particle_flipper(i, round(nucVolume/pow(dx,3)), ranapp);
            SolidD[i] = -nsmooth - 3;
          }
        }
        FDElapsed += time_step;
        time += time_step;
        timer->stamp(TIME_SOLVE);
        if (time >= stoptime)   done = 1;
        if (done || time >= nextoutput) nextoutput = output->compute(time,done);
        timer->stamp(TIME_OUTPUT);
        if (done) break;
      }
      //Compute the normalized mobilities at each site
      //True "mobMax" would be holding the max temp for the entire window

      for(int i = 0; i < nlocal; i++) {
        if(activeFlag[i] < 3) continue;
        MobilityOut[i] = MobilityOut[i]/(mobMax * FDElapsed);
      }
        
    }
    //If FD step is bigger, figure out number of MC steps to do and then update temps.
    //Round to the nearest inter MC steps
    else {
      nMC = round(time_step/dtMC);
      //         if(domain->me == 0) {
      //             fprintf(screen,"Potts time_step %.5e is smaller than solidification time_step %.5e doing %d Potts steps, t_max is %f\n",dtMC, time_step, nMC, tempMaxAll);
      //         }
      //Compute the normalized mobilities for each site
      for(int i = 0; i < nlocal; i++) {
        MobilityOut[i] = compute_mobility(i, ranapp)/mobMax;
      }
    }
  	
  	
    //Do Monte Carlo sweeps
    for (int j = 0; j<nMC; j++) {
      for (int iset = 0; iset < nset_loop; iset++) {
        if (nprocs > 1) {
          timer->stamp();
          if (sectorflag) comm->sector(iset);
          else comm->all();
          timer->stamp(TIME_COMM);
        }

        if (Lmask) boundary_clear_mask(iset);

        timer->stamp();
        // sectors but no colors (could also be no sectors)
        // random selection of sites in iset

        if (sweepflag == RANDOM) {
          site2i = set[iset].site2i;
          nrange = set[iset].nlocal;
          nselect = set[iset].nselect;
          for (i = 0; i < nselect; i++) 
            sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
          (this->*sweep)(nselect,sitelist);
          nattempt += nselect;

          // sectors but no colors, or colors but no sectors
          // ordered sweep over all sites in iset

        } else if (bothflag == 0) {
          //Get rid of nloop because we want this to happen once each time we run it.
          (this->*sweep)(set[iset].nlocal,set[iset].site2i); 
          nattempt += set[iset].nselect;

          // sectors and colors
          // icolor loop is over all colors in a sector
          // jset = set that contains sites of one color in one sector
          // ordered sweep over all sites in jset

        } else {
          for (icolor = 0; icolor < ncolors; icolor++) {
            jset = nsector + iset*ncolors + icolor;
            //Get rid of nloop because we want this to happen once each time we run it.
            (this->*sweep)(set[jset].nlocal,set[jset].site2i);
            nattempt += set[jset].nselect;
          }
        }

        timer->stamp(TIME_SOLVE);

        if (nprocs > 1) {
          if (sectorflag) comm->reverse_sector(iset);
          else comm->all_reverse();
          timer->stamp(TIME_COMM);
        }
      }

      nsweeps++;
      if (dtMC > 1) time += dtMC;
      if (time >= stoptime) done = 1;
      if (done || time >= nextoutput) nextoutput = output->compute(time,done);
      timer->stamp(TIME_OUTPUT);
    }
    //Do a thermal step. Don't do it if dtMC > time_step
    if(dtMC < time_step) {
      app_update(time_step);
      //Check if we're just solidified and should be "relaxed"
      for(int i = 0; i < nlocal; i++) {
            
        if(SolidD[i] < 0 && SolidD[i] > -nsmooth -1)    {
          MobilityOut[i] = 1;
          site_event_rejection(i, ranapp);
          SolidD[i]--;
        }
      }
      //Find the highest temperature in the local array and compare with others
      tempMax = compute_tempMax();
      timer->stamp(TIME_APP);
      MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
      timer->stamp(TIME_COMM);
      //Using the global maximum temperature, compute our smallest timestep
      dtMC = compute_timeMin(tempMaxAll);
      //Also using this value, compute the maximum solid-state Mobility
      mobMax = exp(-Q/(R*tempMaxAll));
      timer->stamp(TIME_APP);
    }
  }
}

/* ----------------------------------------------------------------------
   Function to determine the highest temperature in a region and compute the corresponding
   timestep. Let's compute timestep after and just return highest temperature
   ------------------------------------------------------------------------- */

double AppAdditiveThermal::compute_tempMax() {
  double tempMax = 0;
  //double timeMin = 0;
	
  for(int i = 0; i < nlocal; i++) {
    //Don't include inactive sites in calculation
    if(activeFlag[i] < 2) continue;
    tempMax = MAX(tempMax, T[i]);
    //If max temp is above liquidus, just make it liquidus and return
    if(tempMax > Tl) {
      tempMax = Tl;
      return tempMax;
    }
  }
  return tempMax;
}

/* ----------------------------------------------------------------------
   Given a maximum temperature, compute the minimum MC timestep.
   ------------------------------------------------------------------------- */

double AppAdditiveThermal::compute_timeMin(double tempMax) {
  double dtMC;
  dtMC = pow(dx,2) * Kmc/Ko * exp(Q/(R*tempMax));
  return dtMC;
}
