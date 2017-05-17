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

#include <iostream>
#include <vector>

using std::vector;

namespace RASTER {

namespace pool_shape {

   enum ShapeType {undefined=-1, ellipsoid=0};

   class PoolShape {
      public:
         virtual double distance(const double *XYZ) const = 0;
         PoolShape() {}
         virtual ~PoolShape() {}
   };

}

class point {

public:
	inline point() {p[0]=0.0;p[1]=0.0;p[2]=0.0;}
	inline point(double value) {p[0]=value; p[1]=value; p[2]=value;}
	inline point(double a, double b, double c) { p[0]=a; p[1]=b; p[2]=c; }
	inline point(const double q[3]) { p[0]=q[0]; p[1]=q[1]; p[2]=q[2]; }

	inline point& operator=(const point& rhs) {
      if(this != &rhs){
         p[0]=rhs.p[0];
         p[1]=rhs.p[1];
         p[2]=rhs.p[2];
      }
      return *this;
   }

	inline point(const point& copyMe) {
         p[0]=copyMe.p[0];
         p[1]=copyMe.p[1];
         p[2]=copyMe.p[2];
   }

	inline double squared() const {
		double x=p[X],y=p[Y],z=p[Z];
		return x*x+y*y+z*z;
	}

	inline double operator[](int c) const { return p[c%3]; }

  friend std::ostream& operator<<(std::ostream &os, const point &p)  {
    os << p[X] << ", " << p[Y] << ", " << p[Z] << std::endl;
    return os;
  }

private:
   enum Component {NONE=-1,X=0,Y=1,Z=2};
	double p[3];
};

enum DIR {NONE=-1,X=0,Y=1};

class Pass {
   private:
      DIR _dir;
      double _distance, _speed;
   public:
      Pass(): _dir(DIR::NONE), _distance(-1), _speed(-1) {}
      Pass(DIR d, double distance, double speed): _dir(d), _distance(distance), _speed(speed) {}
      DIR get_dir() const { return _dir; }
      double get_distance() const { return _distance; }
      double get_speed() const { return _speed; }
};

class TransversePass {
   private:
      double _distance, _increment;
   public:
      TransversePass(): _distance(-1), _increment(-1) {}
      TransversePass(double distance, double increment): _distance(distance), _increment(increment) {}
      double get_distance() const { return _distance; }
      double get_increment() const { return _increment; }
};

class RectangularLayer {

   private:
      point _start; 
      double _pool_position[3];
      DIR _dir;
      double _speed;
      double _pass_distance, _overpass, _transverse_pass_distance, _transverse_pass_increment;
      double _pass_distance_traveled, _transverse_pass_distance_traveled;
      int _num_pass;
      bool _serpentine, _on_serpentine;
      // cos(theta), sin(theta)
      double c,s;

   public:
      RectangularLayer()
      : _start(), _pool_position{0,0,0}, _dir(X), _speed(0.0), _pass_distance(0.0), _overpass(0.0), 
        _transverse_pass_distance(0.0), _transverse_pass_increment(0.0), 
        _pass_distance_traveled(0.0), _transverse_pass_distance_traveled(0.0),
        _num_pass(0), _serpentine(false), _on_serpentine(false), 
        c(1.0), s(0.0) {} 

      /*
       * NOTE
       * Layer moves pool up and to the right.
       * Layer initial position is taken as lower left corner.
       */
      RectangularLayer
         (
          const point& start, 
          DIR dir,  
          double speed, 
          double pass_distance, 
          double overpass, 
          double transverse_pass_distance, 
          double transverse_pass_increment,
          bool serpentine=false
         )
      : _start(start), _pool_position{0,0,0}, _dir(dir), _speed(speed), _pass_distance(pass_distance), _overpass(overpass), 
        _transverse_pass_distance(transverse_pass_distance), _transverse_pass_increment(transverse_pass_increment), 
        _pass_distance_traveled(0.0), _transverse_pass_distance_traveled(0.0),
        _num_pass(0), _serpentine(serpentine), _on_serpentine(false),
        c(1.0), s(0.0) 
         {
            // Compute 'cos(theta)' and 'sin(theta)' based upon input pass orientation
            switch(dir){
              case X:
                 c=1.0; s=0.0;
                 break;
              case Y:
                 c=0.0; s=1.0;
                 break;
              default:
                 break;
                 // throw exception here
            }

         } 

      DIR get_dir() const { return _dir; }
      double get_speed() const { return _speed; }
      bool get_serpentine() const {return _serpentine;}
      int get_num_pass() const { return _num_pass;}

      bool move(double dt) {

         bool moved=false;

         double dx=_speed*dt;

         if(_pass_distance_traveled+dx<_pass_distance+_overpass){

            // Update distance traveled
            _pass_distance_traveled+=dx;

            // Increment pool position; pay attention to whether or not we are on 
            //    a serpentine 'pass'
            if(_on_serpentine) _pool_position[0]-=dx; else _pool_position[0]+=dx;

            // Signal that point moved
            moved=true;

         } else if(_transverse_pass_distance_traveled+_transverse_pass_increment<=_transverse_pass_distance){

            // Increment pass counter
            _num_pass+=1;

            // Reset pass distance traveled
            _pass_distance_traveled=0.0;

            // If serpentine, then set direction for serpentine path reversal at odd intervals
            if(_serpentine)  _on_serpentine = (1 ==_num_pass%2);

            // If NOT serpentine, reset  pool_position for travel along pool path
            if(!_serpentine) _pool_position[0]=0.0;
            // If serpentine, we need to subtract the overpass from pool_position, so that
            // we start at the equivalent "_start" everytime
            else if (_on_serpentine)  _pool_position[0]-=_overpass;
            else if (!_on_serpentine) _pool_position[0]+=_overpass;

            // Move in transverse direction
            _transverse_pass_distance_traveled+=_transverse_pass_increment;
            _pool_position[1]=_transverse_pass_distance_traveled;

            // Signal that point moved
            moved=true;
         }
         return moved;
      }

      /*
       * Computes current pool position with respect to 'spparks' 
       * global coordinate system
       *
       */
      point get_position() const {
         // Distance travel with respect to pool coordinates
         double d=_pool_position[0];
         double t=_pool_position[1];
         double dz=0.0;
         
         // Transform above distance traveled to 'spparks' coordinate system
         double dr[]={0.0,0.0,0.0};
         dr[0]= d*c-t*s;
         dr[1]= d*s+t*c;
         dr[2]= dz;

         // Calculate and return new position with respect to 'spparks' coordinate system
         double r[]={0.0,0.0,0.0};
         r[0]=_start[0]+dr[0];
         r[1]=_start[1]+dr[1];
         r[2]=_start[2]+dr[2];
         return point(r);
      }

      /*
       * Compute position of spparks site relative to pool position; Result
       * is expressed with respect to pool local coordinates.  Pool is 
       * always assumed to move along its local x-coordinate axis.
       *
       * xyz: spparks site coordinates triplet (x,y,z)
       * layer_z: z-coordinate of layer managed elsewhere
       *
       */
      point compute_position_relative_to_pool(const double *xyz, double layer_z) const {

         // Pool position with respect to spparks coordinate system
         point pool=get_position();

         // Relative position vector of 'xyz' site coordinate with respect SPPARKS coordinate system
         double x=xyz[0]-pool[0];
         double y=xyz[1]-pool[1];
         double z=xyz[2]-layer_z;

         // Express above position with respect to pool local coordinate system
         double xp= x*c+y*s;
         double yp=-x*s+y*c;
         double zp=z;

         // If on serpentine pass, then negate to properly orient site with respect to 
         //    pool rotated 180 degrees
         if(_on_serpentine){
            xp=-xp;
            yp=-yp;
         }

         return point(xp,yp,zp);
      }
};

class Pattern {
   private:
      int _num_layers;
      vector<int> _layer_ids;
      double _z0, _dz, _z;
      int _current_layer_index=0;
      
   public:
      Pattern() : _num_layers(0), _layer_ids(), _z0(0), _dz(0), _z(0.0), _current_layer_index(0) { } 
      Pattern(const vector<int> layer_ids, double z0, double dz) : 
         _num_layers(layer_ids.size()), _layer_ids(layer_ids), _z0(z0), _dz(dz), _z(z0), _current_layer_index(0) { } 

      double get_layer_z_elevation() const { return _z; }

      int begin() {
         // Increment position to '1' since returning first value;
         _current_layer_index=1;
         // Initialize current elevation
         _z=_z0;
         // return beginning index; 
         return _layer_ids[0];
      }

      int next() {
         int index;
         if(_num_layers==_current_layer_index) {
            // Restart; 
            // Increment position to '1' since returning first value;
            // return beginning index; 
            index=0;
            _current_layer_index=1;

         }
         else {
            index=_current_layer_index;
            _current_layer_index+=1;
         }

         // Increment layer elevation
         _z+=_dz;

         return _layer_ids[index];
      }

};

}

#endif

/*
class RotatedRectangularLayer {
   RotatedRectangularLayer(point<double> start, point<double> end, double angle, double start_offset,  double speed,  double overpass, double orthogonal_increment) {}
};
*/

