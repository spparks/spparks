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

#ifdef APP_CLASS
AppStyle(additive_temperature,AppAdditiveTemperature)

#else

#ifndef SPK_APP_ADDITIVE_TEMPERATURE_H
#define SPK_APP_ADDITIVE_TEMPERATURE_H

#include "app_potts.h"
#include "app_lattice.h"
#include <stdlib.h>

namespace SPPARKS_NS {

class AppAdditiveTemperature : public AppPotts {
 public:
  AppAdditiveTemperature(class SPPARKS *, int, char **);
	virtual void grow_app();
	virtual void init_app();
	virtual void site_event_rejection(int, class RandomPark *);
	virtual double compute_mobility(int, class RandomPark *);
	virtual void nucleation_particle_flipper(int, int,class RandomPark *);
  virtual void input_app(char *, int, char **);
	virtual void app_update(double);
	virtual void mushy_phase(int, class RandomPark *);
	virtual void nucleation_spins(class RandomPark *);
	virtual void nucleation_init();
	virtual void iterate_rejection(double);
	virtual double compute_tempMax();
	virtual double compute_timeMin(double);
	
	
 protected:
// 	int NXeffective;
// 	int NYeffective;
// 	int NZeffective;
// 	double XEllipsePoint;
// 	double YEllipsePoint;
// 	double ZEllipsePoint;
// 	int CurrLayer;
// 	int CurrPass;
	double *MobilityOut;
	double vel;
//     int time_index;
    double x_meltspot, y_meltspot, z_meltspot;
    double minimum_radius;
    double hatch_spacing;
    int end_of_pass_flag; //This variable checks whether the meltspot has completed a total pass length. If so, we ignore the region above/below the YCent depnding on counterclockwise/clockwise
 	
 	//To help improve ease of visualization, lets introduce another integer array. It will be 0 every to begin with
 	//When its layer becomes "active" we'll switch it to 1. This will help visualization and image dumping.
 	int *activeFlag;
	/// parameters for the thermal diffusion eq
	double *T;
	double *SolidD;
// 	double *flux;

	double k_solid;
	double k_powder;
	double time_step;
// 	int diffuse_switch;
// 	double totalTime; //Keep track of time for app_update
// 	int fully_periodic;
	int substrate_height;
	double boundary_temp;
// 	int xPeriod;
// 	int yPeriod;
// 	int zPeriod;
// 	bool initialize_values;         // flag to indicate whether to initialize equilibrium values
	int nlocal_app; //size of cnew, not sure if needed
	
	//Private functions
	void site_event_finitedifference(int);
// 	double boundary_conditioner(int);
// 	double volumetric_gaussian(int);
	double flux_finder(int);
	double specificHeatCalculator(double);
// 	void init_values();
	
	
	double x_dev, y_dev, z_dev; //gaussian heat source size parameters
// 	double mobMax;
	
	//G-Code stuff
	void position_finder_in();
	void path_file();
	void path_file_update();
// 	int gcode_flag;
	  char *path_file_name;
  unsigned line_count;
	double *scan_array;
  double *x_scan_array;
  double *y_scan_array;
  double *z_scan_array;
  double *d_scan_array;
  double *p_scan_array;
//   int FD_steps;
  
  //Let's introduce three new variables to try to make location_finder work better
  double d_residual; //How far have we gone past the last index point
  int path_index; //What index we should start the function looking at
  double wait_time; //How long we should stay in a "wait" block with no temperature input
  
  //We should also define a couple of user input variables to read init_app
  double recoating_time; //How long we should pause the laser between layers
  double short_wait_time; //How long we should wait while the laser is skipping around within a layer
  double flux_prefactor; //The coefficient before the Gaussian heat source
  
	//New parameter stuff
	double density;
	double T_room; //Surrounding temperature for convection/radiation BC's, same as T_inf
	double h; //Convection coefficient
	double eta; //Combined emissivity and stefan boltzmann constant
	
	//New variables for parameterized MC model
	int *nucleationFlags;
	double *nucleationTemps;
	double *nucleationSizes;
	double *neighDist;
	double *specific_heat_temps;
	double *specific_heat_vals;
	int    specific_heat_length;
	double *solid_front_coeffs;
	int    solid_front_length;
	double latent_heat;
	double sizeNorm;
	double sizeSig;
	//New inputs
	double No;
	double Tc;
	double Tsig;
	double dx;
	double Ko; //Scaling pre-factor
	double Kmc; //MC Scaling pre-factor
 	double Q; //Arrhenius exponential factor
 	double Tl; //Liquidus point
 	double Ts; //Solidus point
 	int nsmooth; //How many MC steps to perform after solidification to smooth thigns out
 	const double R = 8.31446261815324; //Define a constant gas constant
 	
 	//Additional variables for parameterized models
	double dtMC; //We should set dt_sweep to this variable at some point
// 	double dtMCMax; //Value at whatever maximum temperature we decide.
// 	double TMCMax;
// 	double MMax; //maximum mobility at the current timestep, set with dtMC
	

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
