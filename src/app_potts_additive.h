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
AppStyle(potts/additive,AppPottsAdditive)

#else

#ifndef SPK_APP_POTTS_ADDITIVE
#define SPK_APP_POTTS_ADDITIVE

#include "app_potts.h"
#include "am_raster.h"
#include <stdlib.h>
#include <map>

using std::map;
using RASTER::Pass;
using RASTER::TransversePass;
using RASTER::RectangularLayer;
using RASTER::Pattern;

namespace SPPARKS_NS {

class AppPottsAdditive : public AppPotts {
 public:
  AppPottsAdditive(class SPPARKS *, int, char **);
  virtual void grow_app();
  virtual void init_app();
  virtual void site_event_rejection(int, RandomPark *);
  void input_app(char *, int , char **);
  double compute_mobility(int, double);
  void app_update(double);
	
 //Remove all of the variables we don't actually need
 protected:
  int NXeffective;
  int NYeffective;
  int NZeffective;
  int CurrLayer;
  int CurrPass;
  
  double *MobilityOut;

 double spot_width;
 double melt_tail_length;
 double tail_HAZ;
 double exp_factor;
 double cap_height;
 double HAZ;
 double melt_depth;
 double cap_HAZ;
 double depth_HAZ;

 private:
   map<int,Pass> passes;
   map<int,TransversePass> transverse_passes;
   map<int,RectangularLayer> rectangular_layers;
   Pattern pattern;
   RectangularLayer active_layer;
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
