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

#ifndef SPK_AM_RASTER_H
#define SPK_AM_RASTER_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>

using std::vector;

namespace RASTER {

namespace pool_shape {

   enum ShapeType {undefined=-1, ellipsoid=0, teardrop=1};

   class PoolShape {
      public:
         virtual double distance(const double *XYZ) const = 0;
         PoolShape() {}
         virtual ~PoolShape() {}
   };

}

class Point {

public:
	inline Point() {p[0]=0.0;p[1]=0.0;p[2]=0.0;}
	inline Point(double value) {p[0]=value; p[1]=value; p[2]=value;}
	inline Point(double a, double b) { p[0]=a; p[1]=b; p[2]=0; }
	inline Point(double a, double b, double c) { p[0]=a; p[1]=b; p[2]=c; }
	inline Point(const double q[3]) { p[0]=q[0]; p[1]=q[1]; p[2]=q[2]; }

	inline Point& operator=(const Point& rhs) {
      //printf("am_rasterh.h; Point& operator=(const Point& p) invoked.\n");
      if(this != &rhs){
         p[0]=rhs.p[0];
         p[1]=rhs.p[1];
         p[2]=rhs.p[2];
      }
      return *this;
   }

	inline Point(const Point& copyMe) {
         //printf("am_rasterh.h; Point(const Layer&p) copy constuctor invoked\n");
         p[0]=copyMe.p[0];
         p[1]=copyMe.p[1];
         p[2]=copyMe.p[2];
   }

	inline double squared() const {
		double x=p[X],y=p[Y],z=p[Z];
		return x*x+y*y+z*z;
	}

	inline double operator[](int c) const { return p[c%3]; }

  friend std::ostream& operator<<(std::ostream &os, const Point &p)  {
    os << p[X] << ", " << p[Y] << ", " << p[Z] << std::endl;
    return os;
  }

private:
   enum Component {NONE=-1,X=0,Y=1,Z=2};
	double p[3];
};

enum DIR {NONE=-1,X=0,Y=1};
/*
 * Starting location for CartesianLayer
 * LL=lower left
 * LR=lower right
 * UR=upper right
 * UL=upper left
 */
enum START {UNDEFINED=-1,LL=0,LR=1,UR=2,UL=3};

/**
 *
 * 2D path; defined in x-y plane by 'start' and 'end' input points
 */
class Path {
   private:
      Point start, end, unit_dir;
      double distance, speed;
   public:

      Path(): start(), end(), distance(0), speed(0) {}

      Path(const Path& p) {
         //printf("am_rasterh.h; Path(const Layer&p) copy constuctor invoked\n");
         start=p.get_start();
         end=p.get_end();
         unit_dir=p.get_unit_dir();
         distance=p.get_distance();
         speed=p.get_speed();
      }

      Path& operator=(const Path& p){
         //printf("am_rasterh.h; Path& operator=(const Path& p) invoked.\n");
         if(&p != this){
            start=p.get_start();
            end=p.get_end();
            unit_dir=p.get_unit_dir();
            distance=p.get_distance();
            speed=p.get_speed();
         }
         return *this;
      }

      Path(const Point& a, const Point& b, double velocity): 
         start(a), end(b), unit_dir(), distance(0.0), speed(velocity) {
    	  Point dr(b[0]-a[0],b[1]-a[1],0);
    	  distance=std::sqrt(dr.squared());
        unit_dir=Point(dr[0]/distance,dr[1]/distance);
      }

      Point get_start() const { return start; }

      Point get_end() const { return end; }

      Point get_unit_dir() const { return unit_dir; }

      double get_distance() const { return distance; }

      double get_speed() const { return speed; }
};

class Pass {
   private:
      DIR dir;
      double speed, hatch_spacing, overhatch;
   public:
      Pass(): dir(DIR::NONE), speed(0), hatch_spacing(0.0), overhatch(0.0) {}
      Pass(DIR d, double vel, double h, double oh): 
         dir(d), speed(vel),  hatch_spacing(h), overhatch(oh) {}
      DIR get_dir() const { return dir; }
      double get_speed() const { return speed; }
      double get_hatch_spacing() const { return hatch_spacing; }
      double get_overhatch() const { return overhatch; }
};

class Layer {
   private:
      int thickness;
      double distance_traveled;
      Point position;
      vector<Path> paths;
      int active_path_ptr;

   public:
      Layer() : thickness(0), distance_traveled(0.0), position(), paths(), active_path_ptr(-1) {}

      Layer(const vector<Path>& _paths, int _thickness) :
         thickness(_thickness), 
         distance_traveled(0.0), 
         position(_paths[0].get_start()),
         paths(_paths),
         active_path_ptr(0) {
            double path_distance=paths[active_path_ptr].get_distance();
            if(path_distance<=0.0){
               throw "Layer constructor; am_raster.h; Invalid path specified; path distance <= 0.0";
            }
      }

      Layer(const Layer& c) : thickness(c.thickness), distance_traveled(c.distance_traveled), 
                              position(c.get_position()), paths(c.paths), 
                              active_path_ptr(c.active_path_ptr) 
      { 
         //printf("Layer(const Layer&c) copy constuctor invoked\n");
         thickness=c.get_thickness();
         distance_traveled=c.distance_traveled;
         position=c.get_position();
         paths=c.paths;
         active_path_ptr=c.active_path_ptr;
      }
      
      Layer& operator=(const Layer& rhs) {

         //printf("Layer& operator=(const Layer& rhs) invoked\n");
         if(this != &rhs){
            thickness=rhs.get_thickness();
            distance_traveled=rhs.distance_traveled;
            position=rhs.get_position();
            paths=rhs.paths;
            active_path_ptr=rhs.active_path_ptr;
         }
         return *this;
      }


      double get_time_on_layer(double dt) const {
         int nsteps=0;
         for(vector<Path>::const_iterator p=paths.begin();p!=paths.end();p++){
            double d=p->get_distance();
            double v=p->get_speed();
            nsteps+=std::ceil(d/(v*dt));
         }
         // -1 is used here to round result down
         return (nsteps-1)*dt;
      }

      int get_thickness() const { return thickness; }

      Point get_unit_dir() const {
         return paths[active_path_ptr].get_unit_dir();
      }

      double get_speed() const {
         return paths[active_path_ptr].get_speed();
      }

      double get_distance_traveled() const {
         return distance_traveled;
      }

   	bool move(double dt) {
         /**
          * Moves pool position along current path.  
          * If pool move is successful, function return 'true'; 
          *    Pool position is computed and relative to SPPARKS coordinate
          *    system
          * Else If path distance has already been traveled returns 'false';
          */
         bool moved=false;
         Point unit_dir=paths[active_path_ptr].get_unit_dir();
         // position increment of motion along path
         double dv=dt*paths[active_path_ptr].get_speed();
         distance_traveled+=dv;

         // Is there more room for motion along path
         if(distance_traveled<paths[active_path_ptr].get_distance()){
            // Signal move
            moved=true;
            // Update position
            Point dp(dv*unit_dir[0],dv*unit_dir[1]);
            // Create new position; copy into temporary;
            Point tmp(position[0]+dp[0],position[1]+dp[1]);
            // Assign new position 
            position=tmp;
         } else {
            // No move on previous path; check and prepare for next path
            active_path_ptr++;
            if(paths.size()!=active_path_ptr){
               // Signal move
               moved=true;
               // position pool to start of next path
               position=paths[active_path_ptr].get_start();
               // initialize distance traveled along path to zero
               distance_traveled=0.0;
            } else {
               // All paths on layer exhausted.
               // Re-initialize position to start of layer;
               //    Re-initializing to start allows for re-use
               //    of this layer in a pattern.
               active_path_ptr=0;
               distance_traveled=0.0;
               position=paths[active_path_ptr].get_start();
               moved=false;
            }
         }
         return moved;
      }

      /*
       * Computes current pool position with respect to 'spparks'
       * global coordinate system
       */
      Point get_position() const { return position; }

      /*
       * See "SPPARKS devel 4" notebook March 9 2020
       *
       * Compute position of spparks site relative to pool position; Result
       * is expressed with respect to pool local coordinates.  Pool is
       * always assumed to move along its local x-coordinate axis.
       *
       * xyz: spparks site coordinates triplet (x,y,z)
       * layer_z: z-coordinate of layer managed elsewhere
       *
       *
       */
      Point compute_position_relative_to_pool(const double *xyz, double layer_z) const {

         // Pool position with respect to spparks coordinate system
         Point pool=get_position();

         // Relative position vector of 'xyz' site coordinate with respect SPPARKS coordinate system
         Point rsx(xyz[0]-pool[0],xyz[1]-pool[1],xyz[2]-layer_z);

         // Project relative site position onto pool local coordinate system
         // d: projection onto path direction
         // t: projection onto transverse direction (right hand rule)
         Point unit_dir=paths[active_path_ptr].get_unit_dir();
         double d=unit_dir[0]*rsx[0]+unit_dir[1]*rsx[1];
         Point dsx(d*unit_dir[0],d*unit_dir[1]);
         Point tsx(rsx[0]-dsx[0],rsx[1]-dsx[1]);
         // dot tsx onto transverse direction to get projection onto 't' axis
         double t=-unit_dir[1]*tsx[0]+unit_dir[0]*tsx[1];

         return Point(d,t,rsx[2]);
      }
};

}

#endif
