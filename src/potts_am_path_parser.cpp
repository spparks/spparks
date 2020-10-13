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

// support class used by several additive manufacturing apps

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
#include "lattice.h"
#include "domain.h"
#include "error.h"
#include "random_mars.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include "am_raster.h"
#include "potts_am_path_parser.h"


using namespace SPPARKS_NS;

using RASTER::DIR;

/* ---------------------------------------------------------------------- */

PottsAmPathParser::PottsAmPathParser(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg), random_park(ranmaster->uniform()),
  passes(), paths(), pattern(), cartesian_layer_meta_data(),
  build_layer_z(std::numeric_limits<double>::quiet_NaN()),num_build_layers(-1), 
  build_layer(0)
{
   // This initializes random_park properly on each processor
   double seed = ranmaster->uniform();
   random_park.reset(seed,domain->me,100);
}

/* ---------------------------------------------------------------------- */

void PottsAmPathParser::print_pool_position(const Point& p)
{
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   {
      if (0==my_rank){
         printf("Pool position: (%6.1f,%6.1f,%6.1f)\n",p[0],p[1],p[2]);
      }
   }
}

/* ---------------------------------------------------------------------- */

void PottsAmPathParser::init_app_am()
{
   delete [] sites;
   delete [] unique;
   sites = new int[1 + maxneigh];
   unique = new int[1 + maxneigh];
   dt_sweep = 1.0/maxneigh;

   // Initialize sites which have value=-1 
   //   * uninitialized sites from stitch file
   //   * if value==-1, initialize to random spin value
   int flag = 0;
   for (int i = 0; i < nlocal; i++){
      if (-1==spin[i]) {
         int ran = random_park.irandom(nspins);
         spin[i]=ran;
      }
      if (spin[i] < 0 || spin[i] > nspins) flag = 1;
   }
   int flagall;
   MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
   if (flagall) error->all(FLERR,"One or more sites have invalid values");

   // Initialize layers
   initialize_layers_am();

}

void PottsAmPathParser::initialize_layers_am()
{
   // Number of layers in script
   int num_pattern_layers=pattern.size();

   // If num_build_layers has not been specified then 
   //    it defaults to number of pattern layers
   if(-1==num_build_layers)
      num_build_layers=pattern.size();

   // Initialize active build layer
   build_layer=0;
   {
      // Compute stop time from given layer/path information
      double _stoptime=0.0;
      for(int l=0;l<num_build_layers;l++){
         int m=l%num_pattern_layers;
         _stoptime+=pattern[m].get_time_on_layer(dt_sweep);
      }
      // Setting stoptime
      stoptime=time+_stoptime;
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      if (0==my_rank){
         // Current build layer
         int m=build_layer%num_pattern_layers;
         Point p=pattern[m].get_position();
         printf("PottsAmPathParser::init_app_am; current time = %8.2f\n", time);
         printf("\tRun time = %8.2f; 'stoptime\' = %8.2f\n", _stoptime, stoptime);
         printf("\tInitial pool position = (%6.1f,%6.1f,%6.1f)\n", p[0],p[1],p[2]);
      }

   }

   // check initialization of 'build_layer_z'
   if(std::isnan(build_layer_z)){
      // NOTE: use 'am_build' command to run simulations where there are
      // lattice sites above build plane; Otherwise, this condition assigns
      // build plane to top surface of domain which is the intended 'stitch'
      // use case.
	   // Input did not include "am_build" command
	   // Assume build_layer_z corresponds to zhi
      build_layer_z=domain->boxzhi;
   }
}

/* ----------------------------------------------------------------------
   manage pool position which is a state variable
------------------------------------------------------------------------- */

bool PottsAmPathParser::app_update_am(double dt)
{
   int num_pattern_layers=pattern.size();
   // Current layer
   int m=build_layer%num_pattern_layers;
   bool moved=true;
   int layer_thickness=pattern[m].get_thickness();
   if(pattern[m].move(dt)){
   } else {
      // Goto next layer
      build_layer++;
      // Have we run out of layers? 
      if(num_build_layers==build_layer){
         int my_rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
         if (0==my_rank){
            printf("PottsAmPathParser::app_update_am; current time = %8.2f\n", time);
            printf("\tPATH END! Time out! 'stoptime\' = %8.2f\n", stoptime);
         }
         moved=false;
         return moved;
      }
      // On a new layer; update z-elevation using 'previous' layer thickness;
      build_layer_z+=layer_thickness;
   }
   return moved;
}

/* ---------------------------------------------------------------------- */

Point PottsAmPathParser::compute_position_relative_to_pool(const double *xyz) const
{
   int num_pattern_layers=pattern.size();
   // Current layer
   int m=build_layer%num_pattern_layers;
   return pattern[m].compute_position_relative_to_pool(xyz,build_layer_z);
}

/* ---------------------------------------------------------------------- */

void PottsAmPathParser::parse_am(int narg, char **arg)
{
   if (strcmp(arg[0],"path") == 0) {
      // printf("%s\n", " PottsAmPathParser::input_app \'am_path\'\n");
      if (narg < 9) error->all(FLERR,"Illegal 'am_path' command; wrong num args.");
      int id=std::atoi(arg[1]);
      double x0,y0,x1,y1;
      double speed;
      if(strcmp(arg[2],"start")==0){
         x0=atof(arg[3]);
         y0=atof(arg[4]);
      } else {error->all(FLERR,"Illegal path command. Expected keyword 'start'");}
      if(strcmp(arg[5],"end")==0){
         x1=atof(arg[6]);
         y1=atof(arg[7]);
      } else {error->all(FLERR,"Illegal path command. Expected keyword 'end'");}
      if(strcmp(arg[8],"speed")==0){
         speed=atof(arg[9]);
      } else {error->all(FLERR,"Illegal path command. Expected keyword 'speed'");}
      paths[id]=Path(Point(x0,y0),Point(x1,y1),speed);
   } else if (strcmp(arg[0],"path_layer") == 0) {
      if (narg < 2) error->all(FLERR,"Illegal 'am_path_layer' command; wrong num args.");
      int id=std::atoi(arg[1]);
      int num_paths,p;
      vector<Path> _paths;
      if(strcmp(arg[2],"num_paths")==0){
         num_paths=std::atoi(arg[3]);
         p=5;
         if(strcmp(arg[4],"path_ids")==0){
            for(int i(0);i<num_paths;i++){
               int pid=std::atoi(arg[p]);
               _paths.push_back(paths[pid]);
               p+=1;
            }
         } else {error->all(FLERR,"Illegal path_layer command. Expected keyword 'path_ids'");}
      } else {error->all(FLERR,"Illegal path_layer command. Expected keyword 'num_paths'");}
      int thickness=-1;
      if(strcmp(arg[p],"thickness")==0){thickness=std::atoi(arg[p+1]); }
      else {error->all(FLERR,"Illegal path_layer command. Expected keyword 'thickness'.");}
      pattern.push_back(Layer(_paths,thickness));
   } else if (strcmp(arg[0],"build") == 0) {
      if (narg != 5) error->all(FLERR,"ERROR 'am build'; wrong num args; expected 5 arguments.");
      double z_start;
      if(strcmp(arg[1],"start")==0){
         z_start=std::atof(arg[2]);
      } else {error->all(FLERR,"Illegal am_build command. Expected keyword 'start'");}
      build_layer_z=z_start;
      if(strcmp(arg[3],"num_layers")==0){
         num_build_layers=std::atoi(arg[4]);
      } else {error->all(FLERR,"Illegal am_build command. Expected keyword 'num_layers'");}
   } else if (strcmp(arg[0],"cartesian_layer") == 0) {
      if (narg < 8) error->all(FLERR,"Illegal 'am cartesian_layer' command; wrong num args; minimum 8 args.");
      add_cartesian_layer(narg,arg);
   } else if (strcmp(arg[0],"pass") == 0) {
      if (narg < 8) error->all(FLERR,"Illegal 'am pass' command; wrong num args; minimum 8 args.");
      add_pass(narg,arg);
   }
}

/* ---------------------------------------------------------------------- */

void PottsAmPathParser::add_pass(int narg, char **arg)
{
   // DEBUG
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if (0==my_rank){
      //printf("PottsAmPathParser::add_pass\n");
   }
   // END DEBUG

   int id=std::atoi(arg[1]);
   DIR d;
   double speed(0.0), hatch(0.0), overhatch(0.0);

   if(strcmp(arg[2],"dir")==0){
      if(strcmp(arg[3],"X")==0) d=DIR::X;
      else if(strcmp(arg[3],"Y")==0) d=DIR::Y;
      else {error->all(FLERR,"Illegal pass 'dir' command. Expected 'X|Y.'"); } 
   } else {error->all(FLERR,"Illegal pass command. Expected 'dir.'");}

   if(strcmp(arg[4],"speed")==0){
      speed=std::atof(arg[5]);
   } else {error->all(FLERR,"Illegal pass command. Expected 'speed.'");}

   if(strcmp(arg[6],"hatch")==0){
      hatch=std::atof(arg[7]);
   } else {error->all(FLERR,"Illegal pass command. Expected 'hatch.'");}
   /*
    * OPTIONAL parameters beyond here
    */
   if(narg>8){
      if (narg!=10) {error->all(FLERR,"Illegal pass command. Expected 8 or 10 arguments only.");}
      if(strcmp(arg[8],"overhatch")==0){
         overhatch=std::atof(arg[9]);
      } else {error->all(FLERR,"ERROR \'pass\' command. Expected \'overhatch' keyword.");}
   } 
   passes[id]=Pass(d,speed,hatch,overhatch);
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PottsAmPathParser::add_cartesian_layer(int narg, char **arg)
{
   // DEBUG
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if (0==my_rank){
      //printf(" PottsAmPathParser::add_cartesian_layer\n");
   }
   // END DEBUG

   int id=std::atoi(arg[1]);
   START s;
   Pass pass;
   int thickness;
   bool serpentine=true;
   double ox(0.0),oy(0.0);
   if(strcmp(arg[2],"start")==0){
      if(strcmp(arg[3],"LL")==0) { s=START::LL; }
      else if(strcmp(arg[3],"LR")==0) s=START::LR;
      else if(strcmp(arg[3],"UR")==0) s=START::UR;
      else if(strcmp(arg[3],"UL")==0) s=START::UL;
      else {error->all(FLERR,"Illegal am_cartesian_layer command  Expected 'LL|LR|UR|UL.'");}
   } else {error->all(FLERR,"Illegal am_cartesian_layer command. Expected 'start.'");}

   if(strcmp(arg[4],"pass_id")==0){
      int id=std::atoi(arg[5]);
      pass=passes[id];
   } else {error->all(FLERR,"Illegal am_cartesian_layer command  Expected 'pass_id'");}

   if(strcmp(arg[6],"thickness")==0){
      thickness=std::atoi(arg[7]);
   } else {error->all(FLERR,"Illegal am_cartesian_layer command  Expected 'thickness.'");}

   /*
    * OPTIONAL parameters beyond here; 
    */
   if(narg>8){
      if(10==narg){
         if(strcmp(arg[8],"serpentine")==0){
            int s=std::atoi(arg[9]);
            // default 'serpentine' is 'true'
            if(0==s) serpentine=false;
         } else {error->all(FLERR,"ERROR am cartesian_layer. 10 arguments.  Expected 'serpentine.'");}
      } else if (11==narg){
         if(strcmp(arg[8],"offset")==0){
            ox=std::atof(arg[9]);
            oy=std::atof(arg[10]);
         } else {error->all(FLERR,"ERROR am cartesian_layer. 11 arguments.  Expected 'offset.'");}
      } else if (13==narg){
         if(strcmp(arg[8],"serpentine")==0){
            int s=std::atoi(arg[9]);
            // default 'serpentine' is 'true'
            if(0==s) serpentine=false;
            if(strcmp(arg[10],"offset")!=0){
               error->all(FLERR,"ERROR am cartesian_layer. 13 arguments.  Expected 'offset.'");
            }
            ox=std::atof(arg[11]);
            oy=std::atof(arg[12]);
         } else if(strcmp(arg[8],"offset")==0){
            ox=std::atof(arg[9]);
            oy=std::atof(arg[10]);
            if(strcmp(arg[11],"serpentine")!=0){
               error->all(FLERR,"ERROR am cartesian_layer. 13 arguments.  Expected 'serpentine.'");
            }
            int s=std::atoi(arg[12]);
            if(0==s) serpentine=false;
         }
      }
   }
   // add cartesian layer meta data
   cartesian_layer_meta_data.push_back(CartesianLayerMetaData(pass,s,ox,oy,thickness,serpentine));
   double xlo=domain->boxxlo;
   double xhi=domain->boxxhi;
   double ylo=domain->boxylo;
   double yhi=domain->boxyhi;
   double x0,y0,x1,y1;
   double speed=pass.get_speed();
   int m=0;
   vector<double>hatch=get_hatch(pass,s,ox,oy);
   vector<Path> layer_paths;
   switch(pass.get_dir()){
   case DIR::X:
      if (START::LL==s || START::UL==s){
         x0=xlo+ox; x1=xhi-ox;
      } else {
         x0=xhi+ox; x1=xlo-ox;
      }
      for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         if(serpentine && (m%2==1)){
            // Flip starting point if 'serpentine' and on 'odd' passes
            layer_paths.push_back(Path(Point(x1,*iter),Point(x0,*iter),speed));
         } else {
            layer_paths.push_back(Path(Point(x0,*iter),Point(x1,*iter),speed));
         }
         m++;
      }
      break;
   case DIR::Y:
      if (START::LL==s || START::LR==s){
         y0=ylo+oy; y1=yhi-oy;
      } else {
         y0=yhi+oy; y1=ylo-oy;
      }
      for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         if(serpentine && (m%2==1)){
            // Flip starting point if 'serpentine' and on 'odd' passes
            layer_paths.push_back(Path(Point(*iter,y1),Point(*iter,y0),speed));
         } else {
            layer_paths.push_back(Path(Point(*iter,y0),Point(*iter,y1),speed));
         }
         m++;
      }
      break;
   }

   // Add this layer
   pattern.push_back(Layer(layer_paths,thickness));

}

vector<double>
PottsAmPathParser::
get_hatch(const Pass& p, START s, double offset_x, double offset_y) const
{
   DIR dir=p.get_dir();
   double hatch_spacing=p.get_hatch_spacing();
   double oh=p.get_overhatch();

   //printf("PottsAmPathParser::get_hatch(...)\n");
   //printf("\thatch spacing=%5.1f",hatch_spacing);
   //printf("\toverhatch=%5.1f\n",oh);

   double xlo=domain->boxxlo;
   double xhi=domain->boxxhi;
   double ylo=domain->boxylo;
   double yhi=domain->boxyhi;
   vector<double> hatch;
   bool start=true;
   if(DIR::X==dir){
      //printf("\tDIR::X\n");
      //printf("\txlo=%8.1f,xhi=%8.1f,ylo=%8.1f,yhi=%8.1f\n",xlo,xhi,ylo,yhi);
      double oy=offset_y;
      if (START::LL==s || START::LR==s){
         double y=ylo+oy;
         while(y<=(yhi+oh)){
            //printf("\tpush_back y=%8.1f\n",y);
            hatch.push_back(y);
            y+=hatch_spacing;
         }
      } else if(START::UL==s || START::UR==s){
         double y=yhi+oy;
         while(y>=(ylo+oh)){
            hatch.push_back(y);
            y-=hatch_spacing;
         }
      } 
   } else if(DIR::Y==dir){
      double ox=offset_x;
      if (START::LL==s || START::UL==s){
         double x=xlo+ox;
         while(x<=(xhi+oh)){
            hatch.push_back(x);
            x+=hatch_spacing;
         }
      } else if(START::LR==s || START::UR==s){
         double x=xhi+ox;
         while(x>=(xlo+oh)){
            hatch.push_back(x);
            x-=hatch_spacing;
         }
      } 
   }

   // Requires c++-14
   // Avoids copies and uses move semantics
   return std::move(hatch);
}

vector<Path> 
PottsAmPathParser::
get_layer_paths(CartesianLayerMetaData& meta) const
{
   Pass pass=meta.get_pass();
   double speed=pass.get_speed();
   bool serpentine=meta.get_serpentine();
   START s=meta.get_start();
   double ox,oy;
   std::tie(ox,oy)=meta.get_offset();
   vector<double>hatch=get_hatch(pass,s,ox,oy);
   vector<Path> layer_paths;
   int m=0;
   double xlo=domain->boxxlo;
   double xhi=domain->boxxhi;
   double ylo=domain->boxylo;
   double yhi=domain->boxyhi;
   double x0,y0,x1,y1;
   switch(pass.get_dir()){
   case DIR::X:
      if (START::LL==s || START::UL==s){
         x0=xlo+ox; x1=xhi-ox;
      } else {
         x0=xhi+ox; x1=xlo-ox;
      }
      for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         if(serpentine && (m%2==1)){
            // Flip starting point if 'serpentine' and on 'odd' passes
            layer_paths.push_back(Path(Point(x1,*iter),Point(x0,*iter),speed));
         } else {
            layer_paths.push_back(Path(Point(x0,*iter),Point(x1,*iter),speed));
         }
         m++;
      }
      break;
   case DIR::Y:
      if (START::LL==s || START::LR==s){
         y0=ylo+oy; y1=yhi-oy;
      } else {
         y0=yhi+oy; y1=ylo-oy;
      }
      for(auto iter=hatch.begin();iter!=hatch.end();iter++){
         if(serpentine && (m%2==1)){
            // Flip starting point if 'serpentine' and on 'odd' passes
            layer_paths.push_back(Path(Point(*iter,y1),Point(*iter,y0),speed));
         } else {
            layer_paths.push_back(Path(Point(*iter,y0),Point(*iter,y1),speed));
         }
         m++;
      }
      break;
   }
   return std::move(layer_paths);
}

vector<PottsAmPathParser::ComputationalVolume>
PottsAmPathParser::
get_cvs(const vector<double>& hatch, const Pass& p, START s, int width_haz) const
{
   DIR dir=p.get_dir();
   int xlo=static_cast<int>(domain->boxxlo);
   int xhi=static_cast<int>(domain->boxxhi);
   int ylo=static_cast<int>(domain->boxylo);
   int yhi=static_cast<int>(domain->boxyhi);
   //printf("PottsAmPathParser::get_cvs\n");
   //printf("\txlo=%8d,xhi=%8d,ylo=%8d,yhi=%8d\n",xlo,xhi,ylo,yhi);
   int cv_x0,cv_x1,cv_y0,cv_y1;
   
   // Make HAZ width just slightly larger than input
   int width=static_cast<int>(1+width_haz/2);
   //printf("\twidth_haz=%8d,width=%8d",width_haz,width);
   vector<ComputationalVolume> cvs;
   for(auto h:hatch){
      //printf("\tDIR::X; Hatch h = %8.1f\n",h);
      if(DIR::X==dir){
         cv_x0=xlo;
         cv_x1=xhi;
         cv_y0=std::floor(h-width);
         cv_y1=std::ceil(h+width);
         if (cv_y0<=ylo)
            cv_y0=ylo;
         if (cv_y1>=yhi)
            cv_y1=yhi;
         //printf("\tcv_x0=%8d,cv_x1=%8d,cv_y0=%8d,cv_y1=%8d\n",cv_x0,cv_x1,cv_y0,cv_y1);
      } else if(DIR::Y==dir){
         cv_y0=ylo;
         cv_y1=yhi;
         cv_x0=std::floor(h-width);
         cv_x1=std::ceil(h+width);
         if (cv_x0<=xlo)
            cv_x0=xlo;
         if (cv_x1>=xhi)
            cv_x1=xhi;
      }
      cvs.push_back(ComputationalVolume(cv_x0,cv_x1,cv_y0,cv_y1));
   }
   // Requires c++-14
   // Avoids copies and uses move semantics
   return std::move(cvs);
}

void PottsAmPathParser::
print_paths
(
   const string& filename, 
   int num_layers, 
   int zstart,
   int melt_depth, 
   int width_haz, 
   int depth_haz
) const
{
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if (0==my_rank){
      FILE* fp = std::fopen(filename.c_str(), "w");
      if(!fp) {
         error->all(FLERR,"potts_am_path_parser.cpp; pathgen file opening failed.");
      }
      const vector<CartesianLayerMetaData>& d=cartesian_layer_meta_data;
      int num_cartesian_layers=d.size();
      int z0,z1=zstart,bottom;
      for(int layer=0;layer<num_layers;layer++){
         int m=layer%num_cartesian_layers;
         CartesianLayerMetaData cartesian_layer=d[m];
         int layer_thickness=cartesian_layer.get_thickness();
         z1+=layer_thickness;
         bottom=z1-layer_thickness;
         if(0==layer) z0=zstart;
         else z0=z1-depth_haz;
         const char* fmt="layer window %8d %8d  haz %8d %8d melt_depth %8d\n";
         //printf(fmt,bottom,z1,z0,z1,melt_depth);
         fprintf(fp,fmt,bottom,z1,z0,z1,melt_depth);
         double ox,oy;
         std::tie(ox,oy)=cartesian_layer.get_offset();
         Pass p=cartesian_layer.get_pass();
         START s=cartesian_layer.get_start();
         vector<double>hatch=get_hatch(p,s,ox,oy);
         vector<ComputationalVolume>cvs=get_cvs(hatch,p,s,width_haz);
         vector<Path> layer_paths=get_layer_paths(cartesian_layer);
         for(int lp=0;lp<cvs.size();lp++){
            int x0,x1,y0,y1;
            std::tie(x0,x1,y0,y1)=cvs[lp];
            Point a=layer_paths[lp].get_start();
            Point b=layer_paths[lp].get_end();
            //printf("path start %8.1f %8.1f end %8.1f %8.1f block %8d %8d %8d %8d \n",
            //      a[0],a[1],b[0],b[1],x0,x1,y0,y1);
            fprintf(fp,"path start %8.1f %8.1f end %8.1f %8.1f block %8d %8d %8d %8d \n",
                  a[0],a[1],b[0],b[1],x0,x1,y0,y1);
         }
      }
      std::fclose(fp);
   }
}
