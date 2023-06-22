
#ifndef SPK_AM_BEZIER_H
#define SPK_AM_BEZIER_H

#include <vector>
#include <array>
#include "kdtree.h"

using std::array;
using std::vector;

class AM_Bezier {

   typedef double value_type;
   typedef std::size_t ordinal_type;
   typedef femanica::kdtree<value_type,ordinal_type> my_tree;

   public:
      AM_Bezier();
      AM_Bezier(const AM_Bezier& q) = default;
      AM_Bezier(AM_Bezier&& q) = default;
      AM_Bezier& operator=(const AM_Bezier& p) = default;
      AM_Bezier& operator=(AM_Bezier&& p) = default;

      AM_Bezier(
                const array<double,5>& x,
                const array<double,5>& y,
                const array<double,5>& z,
                double beta1,
                double beta2,
                double width,
                double depth
            );
      void initialize();

      bool inside_outside(const array<double,3>& xyz, const double tolerance) const;
      bool is_outside_pool(const array<double,3>& xyz) const;
      double distance(const array<double,3>& xyz) const;
      double compute_pool_length() const;
      double get_pool_width() const {return width; }
      double get_pool_length() const {return length; }
      double get_pool_depth() const {return depth; }
      double get_length_to_width_ratio() const {return length/width; }
      double get_length_to_depth_ratio() const {return length/depth; }
      double get_umax_width() const {return umax_width; }
      double get_umax_depth() const {return umax_depth; }
      array<double,3> get_pool_bb_center() const {return bb_center;}
      array<double,3> get_pool_bb_dimensions() const {return bb_dimensions;}
      void print() const;
      void print_rho(double u, double v) const;
      virtual ~AM_Bezier ();
   private:
      // private functions

   private:
      // private variables
      // xt, yt: x and y components of control points for top curve
      // xt: xt are also x component of control points for spine curve
      // zs: z component of control points for spine curve
      // vector<double> xt, yt, zs;
      array<double,5> xt{0.0,0.0,0.0,0.0,0.0};
      array<double,5> yt{0.0,0.0,0.0,0.0,0.0};
      array<double,5> zs{0.0,0.0,0.0,0.0,0.0};
      
      // 4th order polynomial coefficients associated with control points
      array<double,5> cx,cy,cz;

      // 3rd order polynomial coefficients associated with control points
      array<double,4> ry,rz;

      // melt pool width and depth
      double width, umax_width;
      double length;
      double depth, umax_depth;

      array<double,3> bb_center, bb_dimensions;

      // Bezier surface points
      // For use with cpp initial guess
      vector<double> search_points;
      vector<array<double,2>> search_coordinates;
      //shared_ptr<my_tree> tree;
      my_tree tree;
};


#endif

