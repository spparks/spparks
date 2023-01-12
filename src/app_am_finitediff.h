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
AppStyle(am/finitediff,AppAMFiniteDiff)

#else

#ifndef SPK_APP_AM_FINITEDIFF_H
#define SPK_APP_AM_FINITEDIFF_H

#include "app_potts.h"
#include "app_lattice.h"
#include <stdlib.h>

namespace SPPARKS_NS {

class AppAMFiniteDiff : public AppPotts {
 public:
  AppAMFiniteDiff(class SPPARKS *, int, char **);
  virtual void input_app(char *, int, char **);
  virtual void grow_app();
  virtual void init_app();
  virtual void iterate_rejection(double);
  virtual void site_event_rejection(int, class RandomPark *);
  virtual void app_update(double);
	
 protected:
  int done_flag;
  double *MobilityOut;
  double vel;
  double x_meltspot, y_meltspot, z_meltspot;
  double minimum_radius;
  double hatch_spacing;

  // this variable checks whether the meltspot has completed a total
  // pass length. If so, we ignore the region above/below the YCent
  // depnding on counterclockwise/clockwise

  int end_of_pass_flag; 
  
  // To help improve ease of visualization, lets introduce another
  // integer array. It will be 0 every to begin with. When its layer
  // becomes "active" we'll switch it to 1. This will help
  // visualization and image dumping.

  int *activeFlag;

  // parameters for the thermal diffusion eq

  double *T;
  double *SolidD;

  double k_solid;
  double k_powder;
  double time_step;
  int substrate_height;
  double boundary_temp;
  int nlocal_app;             // size of cnew, not sure if needed
	
  double x_dev, y_dev, z_dev;   //  gaussian heat source size parameters

  char *path_file_name;
  unsigned line_count;
  double *scan_array;
  double *x_scan_array;
  double *y_scan_array;
  double *z_scan_array;
  double *d_scan_array;
  double *p_scan_array;
  
  // variables to make location_finder work better

  double d_residual;   // how far have we gone past the last index point
  int path_index;      // what index we should start the function looking at
  double wait_time;    // how long to stay in a "wait" block with no temperature input
  
  // user input variables to read init_app

  double recoating_time;   // how long to pause the laser between layers
  double short_wait_time;  // how long to wait while laser is skips around within a layer
  double flux_prefactor;   // coefficient on the Gaussian heat source
  
  // new parameters

  double density;
  double T_room;    // surrounding temperature for convection/radiation BC's, same as T_inf
  double h;         // convection coefficient
  double eta;       // combined emissivity and stefan boltzmann constant
	
  // new variables for parameterized MC model

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

  // new inputs

  double No;
  double Tc;
  double Tsig;
  double dx;
  double Ko;   // Scaling pre-factor
  double Kmc;  // MC Scaling pre-factor
  double Q;    // Arrhenius exponential factor
  double Tl;   // Liquidus point
  double Ts;   // Solidus point
  int nsmooth; // MC steps to perform after solidification to smooth things out
  const double R = 8.31446261815324; // gas constant
 	
  // additional variables for parameterized models
  
  double dtMC;   // set dt_sweep to this variable at some point

  // private methods

  void path_file();
  void site_event_finitedifference(int);
  double specificHeatCalculator(double);
  void position_finder_in();
  double compute_mobility(int, class RandomPark *);
  double flux_finder(int);
  void path_file_update();
  void mushy_phase(int, class RandomPark *);
  void nucleation_particle_flipper(int, int, class RandomPark *);
  void nucleation_spins(class RandomPark *);
  void nucleation_init();
  double compute_tempMax();
  double compute_timeMin(double);
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
