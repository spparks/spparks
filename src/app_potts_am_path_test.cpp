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
      if (narg != 8) error->all(FLERR,"Illegal 'pathgen' command; wrong num args; should have 8 arguments.");
      if(strcmp(arg[0],"outfile")==0){
         char *filename=arg[1];
         path_filename=std::string(filename);
      } else error->all(FLERR,"Unrecognized command: Expected 'outfile'.");
      if(strcmp(arg[2],"num_layers")==0){
         num_layers=atoi(arg[3]);
      }else error->all(FLERR,"Unrecognized command: Expected 'num_layers'.");
      if(strcmp(arg[4],"melt_depth")==0){
         melt_depth=atoi(arg[5]);
      }else error->all(FLERR,"Unrecognized command: Expected 'melt_depth'.");
      if(strcmp(arg[6],"depth_haz")==0){
         depth_haz=atoi(arg[7]);
      }else error->all(FLERR,"Unrecognized command: Expected 'depth_haz'.");
      // 'pathgen' command must come after all 'am' commands in script
      print_path();
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


void AppPottsAmPathTest::print_path() {

   std::cout << "AppPottsAmPathTest::print_path()" << std::endl;
   std::cout << "\tpathgen outfile = " << path_filename << std::endl;
   std::cout << "\tnum_layers = " << num_layers << std::endl;
   std::cout << "\tmelt_depth = " << melt_depth << std::endl;
   std::cout << "\tdepth_haz = " << depth_haz << std::endl;
   vector<CartesianLayerMetaData> d=get_cartesian_layer_metadata();
   int num_cartesian_layers=d.size();
   for(int layer=0;layer<num_layers;layer++){
         int m=layer%num_cartesian_layers;
         CartesianLayerMetaData cartesian_layer=d[m];
         double ox,oy;
         std::tie(ox,oy)=cartesian_layer.get_offset();
         std::cout << "Layer: " << std::endl;
         std::cout << "\tThickness = " << cartesian_layer.get_thickness() << std::endl;
         //std::cout << "\tPass = " << cartesian_layer.get_pass() << std::endl;
         std::cout << "\tSTART = " << cartesian_layer.get_start() << std::endl;
         std::cout << "\tOffset ox,oy = " << ox << "," << oy << std::endl;
   }
   {
      double xlo=domain->boxxlo;
      double xhi=domain->boxxhi;
      double ylo=domain->boxylo;
      double yhi=domain->boxyhi;
      double x0,y0,x1,y1;
      for(auto layer:d){
         Pass pass=layer.get_pass();
         START s=layer.get_start();
         double ox,oy;
         std::tie(ox,oy)=layer.get_offset();
         double speed=pass.get_speed();
         int m=0;
         vector<double>hatch;
         vector<ComputationalVolume> cvs;
         //tie(hatch,cvs)=get_layer_computational_volumes(pass,s,ox,oy,depth_haz);
         //vector<Path> layer_paths;
         //switch(pass.get_dir()){
         //case DIR::X:
         //   if (START::LL==s || START::UL==s){
         //      x0=xlo+ox; x1=xhi-ox;
         //   } else {
         //      x0=xhi+ox; x1=xlo-ox;
         //   }
         //   for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         //      if(serpentine && (m%2==1)){
         //         // Flip starting point if 'serpentine' and on 'odd' passes
         //         layer_paths.push_back(Path(Point(x1,*iter),Point(x0,*iter),speed));
         //      } else {
         //         layer_paths.push_back(Path(Point(x0,*iter),Point(x1,*iter),speed));
         //      }
         //      m++;
         //   }
         //   break;
         //case DIR::Y:
         //   if (START::LL==s || START::LR==s){
         //      y0=ylo+oy; y1=yhi-oy;
         //   } else {
         //      y0=yhi+oy; y1=ylo-oy;
         //   }
         //   for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         //      if(serpentine && (m%2==1)){
         //         // Flip starting point if 'serpentine' and on 'odd' passes
         //         layer_paths.push_back(Path(Point(*iter,y1),Point(*iter,y0),speed));
         //      } else {
         //         layer_paths.push_back(Path(Point(*iter,y0),Point(*iter,y1),speed));
         //      }
         //      m++;
         //   }
         //   break;
         //}
      }
   }
}

tuple<vector<double>,vector<AppPottsAmPathTest::ComputationalVolume>>
AppPottsAmPathTest::
get_layer_computational_volumes(const Pass& p, START s, int offset_x, int offset_y, int width_haz) const
{
   DIR dir=p.get_dir();
   double hatch_spacing=p.get_hatch_spacing();
   double oh=p.get_overhatch();
   printf("\thatch spacing=%5.1f",hatch_spacing);
   printf("\toverhatch=%5.1f\n",oh);

   double xlo=domain->boxxlo;
   double xhi=domain->boxxhi;
   double ylo=domain->boxylo;
   double yhi=domain->boxyhi;
   vector<double> hatch;
   vector<ComputationalVolume> cvs;
   double cv_x0,cv_x1,cv_y0,cv_y1;
   bool start=true;
   if(DIR::X==dir){
      double oy=offset_y;
      cv_x0=xlo;
      cv_x1=xhi;
      if (START::LL==s || START::LR==s){
         double y=ylo+oy;
         while(y<=(yhi+oh)){
            hatch.push_back(y);
            if(true==start){
               cv_y0=ylo;
               start=false;
            } else {
               cv_y0=y-hatch_spacing;
            }
            y+=hatch_spacing;
            if(y>yhi)
               cv_y1=yhi;
            else
               cv_y1=y;
            // force integer values for block
            cv_y0=std::floor(cv_y0);
            cv_y1=std::ceil(cv_y1);
            cvs.push_back(ComputationalVolume(cv_x0,cv_x1,cv_y0,cv_y1));
         }
      } else if(START::UL==s || START::UR==s){
         double y=yhi+oy;
         while(y>=(ylo+oh)){
            hatch.push_back(y);
            if(true==start){
               cv_y1=yhi;
               start=false;
            } else {
               cv_y1=y+hatch_spacing;
            }
            y-=hatch_spacing;
            if(y<ylo)
               cv_y0=ylo;
            else
               cv_y0=y;
            // force integer values for block
            cv_y0=std::floor(cv_y0);
            cv_y1=std::ceil(cv_y1);
            cvs.push_back(ComputationalVolume(cv_x0,cv_x1,cv_y0,cv_y1));
         }
      } 
   } else if(DIR::Y==dir){
      double ox=offset_x;
      cv_y0=ylo;
      cv_y1=yhi;
      if (START::LL==s || START::UL==s){
         double x=xlo+ox;
         while(x<=(xhi+oh)){
            hatch.push_back(x);
            if(true==start){
               cv_x0=xlo;
               start=false;
            } else {
               cv_x0=x-hatch_spacing;
            }
            x+=hatch_spacing;
            if(x>xhi)
               cv_x1=xhi;
            else
               cv_x1=x;
            // force integer values for block
            cv_x0=std::floor(cv_x0);
            cv_x1=std::ceil(cv_x1);
            cvs.push_back(ComputationalVolume(cv_x0,cv_x1,cv_y0,cv_y1));
         }
      } else if(START::LR==s || START::UR==s){
         double x=xhi+ox;
         while(x>=(xlo+oh)){
            hatch.push_back(x);
            if(true==start){
               cv_x1=xhi;
               start=false;
            } else {
               cv_x1=x+hatch_spacing;
            }
            x-=hatch_spacing;
            if(x<xlo)
               cv_x0=xlo;
            else
               cv_x0=x;
            // force integer values for block
            cv_x0=std::floor(cv_x0);
            cv_x1=std::ceil(cv_x1);
            cvs.push_back(ComputationalVolume(cv_x0,cv_x1,cv_y0,cv_y1));
         }
      } 
   }

   // Requires c++-14
   // Avoids copies and uses move semantics
   return std::make_tuple(std::move(hatch),std::move(cvs));
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
