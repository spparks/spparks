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

#ifndef SPK_POTTS_AM_PATH_PARSER
#define SPK_POTTS_AM_PATH_PARSER

#include <vector>
#include <tuple>
#include <map>
#include <string>
#include "random_park.h"
#include "app_potts.h"
#include "am_raster.h"

using std::map;
using std::vector;
using std::tuple;
using std::string;
using RASTER::Pass;
using RASTER::Path;
using RASTER::Point;
using RASTER::Layer;
using RASTER::START;
using RASTER::CartesianLayerMetaData;

namespace SPPARKS_NS {

class PottsAmPathParser : public AppPotts {

typedef tuple<double,double,double,double> ComputationalVolume;

public:
   PottsAmPathParser(class SPPARKS *, int, char **);
   void print_paths (
      const string& filename, 
      int num_layers, 
      int zstart,
      int melt_depth, 
      int width_haz, 
      int depth_haz) const;

protected:
   void parse_am(int narg, char **arg);
   void init_app_am();
   void initialize_layers_am();
   void print_pool_position(const Point& p);
   bool app_update_am(double dt);
   Point compute_position_relative_to_pool(const double *XYZ) const;

private:
   SPPARKS_NS::RandomPark random_park;
   map<int,Pass> passes;
   map<int,Path> paths;
   vector<Layer> pattern;
   vector<CartesianLayerMetaData> cartesian_layer_meta_data;
   double build_layer_z;
   int num_build_layers, build_layer;
   void add_cartesian_layer(int narg, char **arg);
   void add_pass(int narg, char **arg);
   vector<double>
      get_hatch(const Pass& p, START s, double offset_x, double offset_y) const;
   vector<Path> get_layer_paths(CartesianLayerMetaData& meta) const;
   vector<ComputationalVolume>
      get_cvs(const vector<double>& hatch, const Pass& p, START s, int width_haz) const;
};

}

#endif
