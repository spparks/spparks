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
#include <tuple>
#include <array>
#include <limits>

using std::vector;
using std::tuple;
using std::array;
using std::numeric_limits;

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

private:
   array<double,3> p{0.0,0.0,0.0};

public:
   Point() = default;
   Point(const Point& q) = default;
   Point(Point&& q) = default;
   Point& operator=(const Point& p) = default;
   Point& operator=(Point&& p) = default;

	inline Point(double value) : p{value,value,value} {}
	inline Point(double a, double b) : p{a,b,0.0} {}
	inline Point(double a, double b, double c) : p{a,b,c} {}
	inline Point(const double q[3]) : p{q[0],q[1],q[2]} {}
	inline Point(const array<double,3>& q) : p{q[0],q[1],q[2]} {}
	inline double operator[](int c) const { return p[c%3]; }

	inline double squared() const {
		double x=p[0],y=p[1],z=p[2];
		return x*x+y*y+z*z;
	}

	inline bool operator<=(const Point& r) const {
		return (
            (p[0]<=r[0]) &&
            (p[1]<=r[1]) &&
            (p[2]<=r[2])
            );

	}

	inline bool operator<(const Point& r) const {
		return (
            (p[0]<r[0]) &&
            (p[1]<r[1]) &&
            (p[2]<r[2])
            );
	}


   friend std::ostream& operator<<(std::ostream &os, const Point &p)  {
      os << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
      return os;
   }

};

class rectangular_range {
public:
   typedef double value_type;
	/*
	 * default constructor gives largest possible range
	 */
	rectangular_range()
	:low(-numeric_limits<value_type>::max()),high(numeric_limits<value_type>::max()){}

	rectangular_range(const Point& a, const Point& b)
	: low(a), high(b) {}

	rectangular_range(Point&& a, const Point& b)
	: low(a), high(b) {}

	rectangular_range(const Point& a, Point&& b)
	: low(a), high(b) {}

	rectangular_range(Point&& a, Point&& b)
	: low(a), high(b) {}

	rectangular_range(const Point& c, value_type r)
	: low(c[0]-r,c[1]-r,c[2]-r), high(c[0]+r,c[1]+r,c[2]+r) {}

	rectangular_range& operator=(const rectangular_range&) = default;
	rectangular_range& operator=(      rectangular_range&&) = default;
	rectangular_range(const rectangular_range&) = default;
	rectangular_range(      rectangular_range&&) = default;

	bool contains(const Point& p) const {
		return ((low<=p) && (p<=high));
	}

	bool contains(const rectangular_range&r) const {
		return ((low<=r.get_low()) && (r.get_high() <= high));
	}

	bool intersects(const rectangular_range&r) const {
		return (!(high<r.get_low()) && !(r.get_high()<low));
	}

	Point get_low() const { return low; }
	Point get_high() const { return high; }

	friend std::ostream& operator<<(std::ostream &os, const rectangular_range& r){
		os << r.get_low() << r.get_high();
		return os;
	}

private:
	Point low, high;
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
      Point start={0.0,0.0,0.0}, end={0.0,0.0,0.0}, unit_dir={1.0,0.0,0.0};
      double distance=0.0, speed=0.0;

   public:
      Path(const Point& a, const Point& b, double velocity): 
        start(a), end(b), speed(velocity) {
    	  Point dr(b[0]-a[0],b[1]-a[1],0);
    	  distance=std::sqrt(dr.squared());
        unit_dir=Point(dr[0]/distance,dr[1]/distance);
      }
      Path() = default;
      Path(const Path& q) = default;
      Path(Path&& q) = default;
      Path& operator=(const Path& p) = default;
      Path& operator=(Path&& p) = default;

   public:
      Point get_start() const { return start; }

      Point get_end() const { return end; }

      Point get_unit_dir() const { return unit_dir; }

      double get_distance() const { return distance; }

      double get_speed() const { return speed; }

   friend std::ostream& operator<<(std::ostream &os, const Path &p)  {
      os << "Start: "  << p.get_start() << "\n"
         << "End: " << p.end  << "\n"
         << "Unit Dir: " << p.unit_dir << "\n"
         << "Distance = " << p.distance << "\n"
         << "Speed = " << p.speed << std::endl;
      return os;
   }
};

class Pass {
   private:
      DIR dir=DIR::NONE;
      double speed=0.0, hatch_spacing=0.0, overhatch=0.0;
   public:
      Pass() = default;
      Pass(const Pass& p) = default;
      Pass(Pass&& p) = default;
      Pass& operator=(const Pass& p) = default;
      Pass& operator=(Pass&& p) = default;
      Pass(DIR d, double vel, double h, double oh): 
         dir(d), speed(vel),  hatch_spacing(h), overhatch(oh) {}
      DIR get_dir() const { return dir; }
      double get_speed() const { return speed; }
      double get_hatch_spacing() const { return hatch_spacing; }
      double get_overhatch() const { return overhatch; }
};

class CartesianLayerMetaData {
   private:
      Pass p=Pass();
      START s=START::UNDEFINED;
      double offset_x=0.0, offset_y=0.0;
      int thickness=0;
      bool serpentine=true;
   public:
      CartesianLayerMetaData() = default;
      CartesianLayerMetaData(const CartesianLayerMetaData& m) = default;
      CartesianLayerMetaData(CartesianLayerMetaData&& m) = default;
      CartesianLayerMetaData& operator=(const CartesianLayerMetaData& m) = default;
      CartesianLayerMetaData& operator=(CartesianLayerMetaData&& m) = default;
      CartesianLayerMetaData(
            const Pass q,
            const START start,
            double ox,
            double oy,
            int t,
            bool s
            ):
         p(q), s(start), offset_x(ox), offset_y(oy), thickness(t), serpentine(s) {}

      Pass get_pass() const { return p; }
      START get_start() const { return s; }
      int get_thickness() const { return thickness; }
      bool get_serpentine() const { return serpentine; }
      tuple<double,double> get_offset() const {return std::make_tuple(offset_x,offset_y);}
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
