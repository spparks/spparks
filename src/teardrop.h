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

#ifndef SPK_TEARDROP_H
#define SPK_TEARDROP_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <limits>
#include <stdio.h>
#include <tuple>
#include "pool_shape.h"

namespace weld {

namespace pool_shape {

   using std::vector;

   double dot2d(const double* y1, const double *y2) { return y1[0]*y2[0]+y1[1]*y2[1]; }
   double dot3d(const double* y1, const double *y2) { return y1[0]*y2[0]+y1[1]*y2[1]+y1[2]*y2[2]; }

   /* ---------------------------------------------------------------------- */

   vector<double> evaluate_bernstein_polynomials(int n, double u) {
      vector<double> b(n+1); b[0]=1.0; 
      double u1=1.0-u;
      for(int j=1;j<=n; j++){
        double bj=0.0; 
        for(int k=0;k<j;k++){
           double t=b[k];
           b[k]=bj+u1*t;
           bj=u*t;
        }
        b[j]=bj;
      }

      return b;
   }

   /* ---------------------------------------------------------------------- */
   
   void evaluate_bernstein_polynomials(double u, vector<double>& b) {
      // n is polynomial order and is implied by the length of 'b'
      int n=b.size()-1;
      b[0]=1.0;
      double u1=1.0-u;
      for(int j=1;j<=n; j++){
        double bj=0.0; 
        for(int k=0;k<j;k++){
           double t=b[k];
           b[k]=bj+u1*t;
           bj=u*t;
        }
        b[j]=bj;
      }
   }

   /* ---------------------------------------------------------------------- */

   class BezierCurve {

      private:
         int n,dim;
         vector< vector<double> > P, dP, ddP;

      public:
         BezierCurve(const vector< vector<double> >& _control_points): 
            n(_control_points.size()-1), dim(_control_points[0].size()), 
            P(_control_points), dP(_control_points.size()-1), ddP(_control_points.size()-2) {
            
            // Using input control points, updates/computes 'dP', and 'ddP'
            update_control_point_deltas();
         }

         void update_control_point_deltas(){
            // Each control point is a vector of dim=2 or 3.  Assume that all control points have same dimension.
            {
               // For purposes of calculating first derivative, compute differences of control points.
               vector<double> p, pp1, _dp(dim);
               for(int i=0;i<=n-1;i++){
                  pp1=P[i+1];
                  p=P[i];
                  for(int d=0;d<dim;d++)
                     _dp[d]=pp1[d]-p[d];
                  dP[i]=_dp;
               }
            }
            {
               // For purposes of calculating second derivative, compute differences of differences of control points.
               vector<double> dp, dpp1, _ddp(dim);
               for(int i=0;i<=n-2;i++){
                  dpp1=dP[i+1];
                  dp=dP[i];
                  for(int d=0;d<dim;d++)
                     _ddp[d]=dpp1[d]-dp[d];
                  ddP[i]=_ddp;
               }
            }

         }

         void move_curve(double dy) {
            // This function CHANGES 'this' curve.
            // Moves curve in 'y-direction' by 'dy'.
            // NOTE that 'delta' control points do not have to be 
            //    since a move has no effect on them.
            int N=n;
            for(int i=0;i<=N;i++){
               vector<double> pi=P[i];
               // Only move 'y-coordinates'
               pi[1]+=dy;
               P[i]=pi;
            }
         }

         void scale_curve(double alpha) {
            // This function CHANGES 'this' curve.
            // Scales curve by multiplying control points with 
            //    input 'alpha'
            int N=n;
            for(int i=0;i<=N;i++){
               vector<double> pi=P[i];
               for(int j=0;j<dim;j++)
                  pi[j]*=alpha;
                  P[i]=pi;
            }
            update_control_point_deltas();
         }

         void compute_curve(double u, double *curve) const {
            int N=n;
            double *c=curve; 
            for(int j=0;j<dim;j++)
               c[j]=0.0;
            vector<double> b=evaluate_bernstein_polynomials(N,u);
            for(int i=0;i<=N;i++){
               vector<double> pi=P[i];
               double bi=b[i];
               for(int j=0;j<dim;j++)
                  c[j]+=pi[j]*bi;
            }
         }

         void compute_dcurve(double u, double *dcurve) const {
            // Use Bernstein polynomials of 1 order less than curve 
            //    used to model teardrop shape
            int N=n-1;
            double *dc=dcurve; 
            for(int j=0;j<dim;j++)
               dc[j]=0.0;
            vector<double> b=evaluate_bernstein_polynomials(N,u);
            // Use 'dP' for control points when computing derivative of curve
            // Also note multiplication of Bernstein polynomials by 'n'
            for(int i=0;i<=N;i++){
               vector<double> pi=dP[i];
               double bi=n*b[i];
               for(int j=0;j<dim;j++)
                  dc[j]+=pi[j]*bi;
            }
         }

         void compute_ddcurve(double u, double *ddcurve) const {
            // Use Bernstein polynomials of 2 orders less than curve 
            //    used to model teardrop shape
            int N=n-2;
            double *ddc=ddcurve; 
            for(int j=0;j<dim;j++)
               ddc[j]=0.0;
            vector<double> b=evaluate_bernstein_polynomials(N,u);
            // Use 'ddP' for control points when computing 2nd derivative of curve
            // Also note multiplication of Bernstein polynomials by '(n)*(n-1)'
            for(int i=0;i<=N;i++){
               vector<double> pi=ddP[i];
               double bi=(n)*(n-1)*b[i];
               for(int j=0;j<dim;j++)
                  ddc[j]+=pi[j]*bi;
            }
         }

         int curve_spatial_dimension() const { return dim; }

   };

   /* ---------------------------------------------------------------------- */

   class TeardropCurve {

      private:
         vector< vector<double> > P;
         const int ai=0, bi=1, ci=2, di=3, ei=4;

      public:

         void compute_curve(double u, double *curve) const {
            double *c=curve; 

            const double u2=u*u;
            const double u3=u2*u;
            const double u4=u2*u2;

            const double um1=(1-u);
            const double um1_2=um1*um1;
            const double um1_3=um1_2*um1;
            const double um1_4=um1_2*um1_2;

            const double um1_3u=um1_3*u;
            const double um1_2u2=um1_2*u2;
            const double um1u3=um1*u3;

            const double bx=P[bi][0];
            const double cx=P[ci][0];
            const double dx=P[di][0];

            const double ay=P[ai][1];
            const double cy=P[ci][1];
            const double dy=P[di][1];

            // X-component of curve
            c[0]=4*bx*um1_3u+6*cx*um1_2u2+4*dx*um1u3;

            // Y-component of curve
            c[1]=ay*um1_4+4*ay*um1_3u+6*cy*um1_2u2+4*dy*um1u3+dy*u4;
         }

         void compute_dcurve(double u, double *dcurve) const {
            double *dc=dcurve; 

            const double u2=u*u;
            const double u3=u2*u;

            const double um1=(1-u);
            const double um1_2=um1*um1;
            const double um1_3=um1_2*um1;

            const double um1_2u=um1_2*u;
            const double um1u2=um1*u2;

            const double bx=P[bi][0];
            const double cx=P[ci][0];
            const double dx=P[di][0];

            const double ay=P[ai][1];
            const double cy=P[ci][1];
            const double dy=P[di][1];

            // X-component of derivative of the curve
            dc[0]=4*bx*um1_3+12*(cx-bx)*um1_2u+12*(dx-cx)*um1u2-4*dx*u3;

            // Y-component of derivative of the curve
            dc[1]=12*(cy-ay)*um1_2u+12*(dy-cy)*um1u2;
         }

         void compute_ddcurve(double u, double *ddcurve) const {
            double *ddc=ddcurve; 

            const double u2=u*u;

            const double um1=(1-u);
            const double um1_2=um1*um1;

            const double um1u=um1*u;

            const double bx=P[bi][0];
            const double cx=P[ci][0];
            const double dx=P[di][0];

            const double ay=P[ai][1];
            const double cy=P[ci][1];
            const double dy=P[di][1];

            // X-component of 2nd derivative of the curve
            ddc[0]=(12*cx-24*bx)*um1_2+(24*bx-48*cx+24*dx)*um1u+(12*cx-24*dx)*u2;

            // Y-component of 2nd derivative of the curve
            ddc[1]=12*(cy-ay)*um1_2+(24*ay-48*cy+24*dy)*um1u+12*(cy-dy)*u2;
         }

         void compute_curve(double u, double *curve, double *dcurve) const {
            double *c=curve, *dc=dcurve; 

            const double u2=u*u;
            const double u3=u2*u;
            const double u4=u2*u2;

            const double um1=(1-u);
            const double um1_2=um1*um1;
            const double um1_3=um1_2*um1;
            const double um1_4=um1_2*um1_2;

            const double um1_3u=um1_3*u;
            const double um1_2u2=um1_2*u2;
            const double um1u3=um1*u3;
            const double um1_2u=um1_2*u;
            const double um1u2=um1*u2;

            const double bx=P[bi][0];
            const double cx=P[ci][0];
            const double dx=P[di][0];

            const double ay=P[ai][1];
            const double cy=P[ci][1];
            const double dy=P[di][1];

            // X-component of curve
            c[0]=4*bx*um1_3u+6*cx*um1_2u2+4*dx*um1u3;
            // X-component of derivative of the curve
            dc[0]=4*bx*um1_3+12*(cx-bx)*um1_2u+12*(dx-cx)*um1u2-4*dx*u3;

            // Y-component of curve
            c[1]=ay*um1_4+4*ay*um1_3u+6*cy*um1_2u2+4*dy*um1u3+dy*u4;
            // Y-component of derivative of the curve
            dc[1]=12*(cy-ay)*um1_2u+12*(dy-cy)*um1u2;
         }

         void compute_curve(double u, double *curve, double *dcurve, double *ddcurve) const {
            double *c=curve, *dc=dcurve, *ddc=ddcurve; 

            const double u2=u*u;
            const double u3=u2*u;
            const double u4=u2*u2;

            const double um1=(1-u);
            const double um1_2=um1*um1;
            const double um1_3=um1_2*um1;
            const double um1_4=um1_2*um1_2;

            const double um1_3u=um1_3*u;
            const double um1_2u2=um1_2*u2;
            const double um1u3=um1*u3;
            const double um1_2u=um1_2*u;
            const double um1u2=um1*u2;
            const double um1u=um1*u;

            const double bx=P[bi][0];
            const double cx=P[ci][0];
            const double dx=P[di][0];

            const double ay=P[ai][1];
            const double cy=P[ci][1];
            const double dy=P[di][1];

            // X-component of curve
            c[0]=4*bx*um1_3u+6*cx*um1_2u2+4*dx*um1u3;
            // X-component of derivative of the curve
            dc[0]=4*bx*um1_3+12*(cx-bx)*um1_2u+12*(dx-cx)*um1u2-4*dx*u3;
            // X-component of 2nd derivative of the curve
            ddc[0]=(12*cx-24*bx)*um1_2+(24*bx-48*cx+24*dx)*um1u+(12*cx-24*dx)*u2;

            // Y-component of curve
            c[1]=ay*um1_4+4*ay*um1_3u+6*cy*um1_2u2+4*dy*um1u3+dy*u4;
            // Y-component of derivative of the curve
            dc[1]=12*(cy-ay)*um1_2u+12*(dy-cy)*um1u2;
            // Y-component of 2nd derivative of the curve
            ddc[1]=12*(cy-ay)*um1_2+(24*ay-48*cy+24*dy)*um1u+12*(cy-dy)*u2;
         }

         double get_pool_length() const {
            double ymin,ymax;
            // Length of pool is controlled by y-component 
            //   of first and last control point
            ymin=P[0][1];
            ymax=P[4][1];
            double pool_length=ymax-ymin;
            return pool_length;
         }

         double get_parametric_coordinate_at_maximum_pool_width() const {
            const double u0=0.5;
            double u=u0;
            int iter=0;
            const double tolerance=1.0e-15;
            const int max_iter=15;
            const int dim=2;
            double c[dim],dc[dim],ddc[dim];
            compute_curve(u,c,dc,ddc);
            double f=dc[0];
            double j=ddc[0];
            double resid=std::fabs(f);
            while (resid>tolerance && iter<=max_iter){
               // Compute increment and update newton iteration
               // u+du = inv(Jacobian) * residual
               u+=-f/j;
               compute_curve(u,c,dc,ddc);
               f=dc[0];
               j=ddc[0];
               resid=std::fabs(f);
               iter+=1;
            }
            return u;
         }
   
         double get_pool_width() const {
            // Compute width of pool at parametric coordinate where maximum
            //    pool width occurs; SEE Comments in function which computes umax;
            double umax=get_parametric_coordinate_at_maximum_pool_width();
            const int dim=2;
            double c[dim];
            compute_curve(umax,c);
            double xmax=std::fabs(c[0]);
            double xmin=-xmax;
            double pool_width=xmax-xmin;
            return pool_width;
         }

         double get_pool_position() const {
            // Pool position is defined average of geometric
            //    centroid along y-axis
            double ymin,ymax;
            ymin=P[0][1];
            ymax=P[4][1];
            double pool_position=0.5*(ymax+ymin);
            return pool_position;
         }


         vector< vector<double> > get_copy_of_control_points() const { return P; } 

         ~TeardropCurve() {}
         TeardropCurve(const vector< vector<double> >& control_points, double pool_width, bool normalize_pool_position=false) :
            P(control_points) {
            {
               // Normalize set of control points shape of pool; set 'pool width' ;
               // Pool length will automatically follow according to shape of control points.
               double W0=get_pool_width();
               // Scale pool size to pool_width
               double alpha=pool_width/W0;
               scale_curve(alpha);

               if(normalize_pool_position) {
                  // Positions pool at y=0.0
                  double dy=0.0;
                  double yp=get_pool_position();
                  dy-=yp;
                  move_curve(dy);
               }
            }
         }

      private:

         void scale_curve(double alpha) {
            // This function CHANGES 'this' curve.
            // Scales curve by multiplying control points with 
            //    input 'alpha'
            const int N=4;
            const int dim=2;
            for(int i=0;i<=N;i++){
               vector<double> pi=P[i];
               for(int j=0;j<dim;j++)
                  pi[j]*=alpha;
                  P[i]=pi;
            }
         }

         void move_curve(double dy) {
            // This function CHANGES 'this' curve.
            // Moves curve in 'y-direction' by 'dy'.
            // NOTE that 'delta' control points do not require an update
            //    since a move has no effect on them.
            const int N=4;
            for(int i=0;i<=N;i++){
               vector<double> pi=P[i];
               // Only move 'y-coordinates'
               pi[1]+=dy;
               P[i]=pi;
            }
         }
   
   };

   /* ---------------------------------------------------------------------- */

   class Teardrop2D : public TeardropCurve {

      private:

      public:
         ~Teardrop2D() {}

         Teardrop2D(const vector< vector<double> >& _teardrop_control_points, double pool_width, bool normalize_pool_position) : 
            TeardropCurve(_teardrop_control_points,pool_width,normalize_pool_position) { }


         double predict_cpp_parametric_coordinate(double x, double y, double du) const {

            const double tolerance=1.0e-15;
            double um=get_parametric_coordinate_at_maximum_pool_width();
            double ym;
            {
               // Compute y-coordinate at maximum pool width
               double rho[]={0,0};
               compute_curve(um,rho);
               ym=rho[1];
            }
            double s;
            double umin, umax, minlim, maxlim;
            // Identify quadrant of input spatial point (x,y)
            // Starting value for bounds on parametric coordinate
            minlim=0.0;
            maxlim=1.0;
            if (y>ym){
               umin=um;
               umax=1.0;
               // If point is on weld axis
               if(x<tolerance) return umax;
            }
            else {
               umin=0;
               umax=um;
               // If point is on weld axis
               if(x<tolerance) return umin;
            }

            {
               // Bound parametic coordinate and compute (umin, umax) of bounds
               double un=0.5*(umin+umax);
               while (umax-umin > du){
                  double rho[]={0,0};
                  double rho_1[]={0,0};
                  compute_curve(un,rho,rho_1);
                  double dr[]={x-rho[0],y-rho[1]};
                  double dt=dot2d(dr,rho_1);
                  if (dt>0){
                      // Point is north of un
                      umin=un;
                  }
                  else{
                      // Point is south of un
                      umax=un;
                  }
                  un=0.5*(umin+umax);
               }
            }

            {
               // Compute points on curve at c0(umin), c1(umax)
               double c0[]={0,0};
               double c1[]={0,0};
               compute_curve(umin,c0);
               compute_curve(umax,c1);
               double delta[]={c0[0]-c1[0],c0[1]-c1[1]};
               double dr[]={x-c0[0],y-c0[1]};
               double t=dot2d(delta,dr)/dot2d(delta,delta);
               if(0>t || t >1){
                  // Override and limit value of 't'
                  if(t>1 && umax<maxlim){
                     t=1;
                  } else {
                     if(t<0 && umin>minlim){
                        t=0;
                     } else{ t=0.5; }
                  }
               }
               s=(1-t)*umin+t*umax;
            }
            // Predicted parametric coordinate for closest point projection
            return s;
         }

         std::tuple<double,double> closest_point(double x, double y) const {
            // Compute parametric coordinate defining closest point;  
            // Using parametric coordinate for closest point on curve, 
            // function also return distance to curve.
            //
            // Initial guess must be decent otherwise cpp fails
            double du=0.025;
            double u0=predict_cpp_parametric_coordinate(x,y,du);

            // Initial estimate
            double u=u0;
            double rho[]={0,0};
            double rho_1[]={0,0};
            double rho_11[]={0,0};
            compute_curve(u,rho,rho_1,rho_11);
            double dr[]={x-rho[0],y-rho[1]};
            {
               // Is this point an interior point?  
               // 'rho_1' is tangent vector
               // Define normal 'n' to curve using right hand rule so that it points 
               //    towards interior of pool; 
               double n[]={-rho_1[1],rho_1[0]};
               // Minus sign because of the way 'dr' was defined; dr is 
               //   currently vector pointing from 'curve' to input point 'x,y' 
               //   which is exactly the opposite of what is needed in 
               //   the following dot product. 
               double dot=-(n[0]*dr[0]+n[1]*dr[1]);
               // If above dot product is negative then (x,y) defines 
               // a point inside the pool.
               if (dot<=0) {
                  return std::make_tuple(0.0,-1.0);
               }
            }
            // compute_ddcurve(u,rho_11);
            // Jacobian
            double j11=dot2d(rho_11,dr)-dot2d(rho_1,rho_1);
            //std::cout << "dr = " << std::scientific << std::setw(15) << dr[0] << ", " << dr[1] << std::endl;
            //std::cout << "rho_1 = " << std::scientific << std::setw(15) << rho_1[0] << ", " << rho_1[1] << std::endl;

            // Compute residual
            double r1=dot2d(rho_1,dr);
            double resid=std::fabs(r1);

            int iter=0;
            //std::cout << "Iteration # Solution Value  Residual" << std::endl;
            //std::cout << "----------- --------------  --------" << std::endl;
            //std::cout << std::setw(10) << iter;
            //std::cout << std::scientific << std::setw(15) << u;
            //std::cout << std::scientific << std::setw(15) << resid << std::endl;
            double tolerance=1.0e-15;
            int max_iter=15;
            while (resid>tolerance && iter<=max_iter){

               // Compute increment and update newton iteration
               // u+du = inv(Jacobian) * residual
               u+=-r1/j11;

               // Update residual in preparation for next iteration 
               compute_curve(u,rho,rho_1,rho_11);
               //compute_curve(u,rho);
               dr[0]=x-rho[0];
               dr[1]=y-rho[1];
               // compute_dcurve(u,rho_1);
               r1=dot2d(rho_1,dr);
               resid=std::fabs(r1);

               // In preparation of next Jacobian calculation
               // compute_ddcurve(u,rho_11);
               // Jacobian
               j11=dot2d(rho_11,dr)-dot2d(rho_1,rho_1);

               // Go to next iteration
               iter+=1;
               //std::cout << std::setw(10) << iter;
               //std::cout << std::scientific << std::setw(15) << u;
               //std::cout << std::scientific << std::setw(15) << resid << std::endl;
            }

            // Identify points which are inside pool by setting distance=-1;
            // 'rho_1' is tangent vector
            // Define normal to curve using right hand rule so that it points 
            //    towards interior of pool; 
            double n[]={-rho_1[1],rho_1[0]};
            double distance=std::sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
            // Minus sign because of the way 'dr' was defined; dr is 
            //   currently vector pointing from 'curve' to input point 'x,y' 
            //   which is exactly the opposite of what is needed in 
            //   the following dot product.
            double dot=-(n[0]*dr[0]+n[1]*dr[1]);
            if (dot<=0) distance=-1;
            return std::make_tuple(u,distance);

         }

   };

   /* ---------------------------------------------------------------------- */

   class Teardrop : public TeardropCurve, public PoolShape {

      private:
         const double _h, _alpha, _beta;

      public:
         Teardrop(
                  double plate_thickness, 
                  double pool_width,
                  double alpha, 
                  double beta,
                  const vector< vector<double> >& teardrop_control_points
                  ) : TeardropCurve(teardrop_control_points,pool_width,true), PoolShape(), 
            _h(plate_thickness), _alpha(alpha), _beta(beta)
         {}

         virtual ~Teardrop(){}

         virtual double distance(const double *XYZ) const {
            /**
             * For points 'xyz' outside of pool, 
             * computes closest point projection to 
             * 'Teardrop' surface; Using this 
             * projection the distance to the surface 
             * is returned.
             */

            // Model uses symmetry about y-z plane; 
            // Model assumes and requires x-coordinates are positive.
            // Because of symmetry, distance(x)=distance(-x)
            double xyz[]={std::fabs(*XYZ), *(XYZ+1), *(XYZ+2)};
            // Predict parameteric coordinates for closest point project
            double u0, v0, distance;
            {
               // Compute initial guess (u0,v0) at parametric coordinates
               double x=*xyz;
               double y=*(xyz+1);
               double z=*(xyz+2);
               double z2=z*z;
               double h2=_h*_h;
               double c=(h2+2*_h*z*(-1+_alpha)*(-1+_beta)+z2*(-1+_alpha)*(-1+2*_beta))/h2;
               double W=this->get_pool_width();
               vector< vector<double> > cp=this->get_copy_of_control_points();
               bool normalize_pool_position=false;
               Teardrop2D t2d(cp,c*W,normalize_pool_position);
               std::tie(u0,distance)=t2d.closest_point(x,y);
               if(distance<0.0) return -1.0;

               // SHORT CIRCUIT for 2D case
               // return distance;
               // END SHORT CIRCUIT
               v0=-z/_h;
            }
            double U0[]={u0,v0};
            double U[]={0,0};
            compute_closest_point(xyz,U0,U);

            // Compute distance 'd' here
            double u=*U;
            double v=*(U+1);
            double rho[]={0,0,0};
            if(v<0){
               // If v < 0, then closest point is above pool surface (outside of domain)
               //         In that case, assume closest point is at pool surface --> v=0
               u=u0; v=0;
            }
            compute_rho(u,v,rho);
            double dr[]={xyz[0]-rho[0],xyz[1]-rho[1],xyz[2]-rho[2]};
            double d=sqrt(dot3d(dr,dr));
            return d;
         }

      private:
         void compute_closest_point(const double *xyz, const double *u0, double *U) const {
            /*
             * Computes closest point 'U' on surface for given input point 'xyz' outside pool.
             *
             * Inputs
             * ------
             * All inputs must be pre-allocated to length described below.
             *
             * xyz: 3 component vector for point; seek to find pair (u,v) on pool surface
             *      which is closest point projection of this input.
             *
             * u0: 2 component vector storing initial guess at closest point projection.
             *
             * U: 2 component vector solution givin parametric coordinates of surface.
             *    
             * Outputs
             * -------
             *  Computed closest point U in parametric coordinates of surface.
             * 
             */

            double u=*u0;
            double v=*(u0+1);
            {
               double rho[]={0,0,0};
               compute_rho(u,v,rho);
               double dr[]={xyz[0]-rho[0],xyz[1]-rho[1],xyz[2]-rho[2]};
               double rho_1[]={0,0,0};
               compute_rho_1(u,v,rho_1);
               double rho_2[]={0,0,0};
               compute_rho_2(u,v,rho_2);
               double rho_11[]={0,0,0};
               compute_rho_11(u,v,rho_11);
               double rho_12[]={0,0,0};
               compute_rho_12(u,v,rho_12);
               double rho_21[]={0,0,0};
               compute_rho_21(u,v,rho_21);
               double rho_22[]={0,0,0};
               compute_rho_22(u,v,rho_22);

               // Compute residual
               double r1=dot3d(rho_1,dr);
               double r2=dot3d(rho_2,dr);
               double resid=sqrt(r1*r1+r2*r2);

               double tolerance=1.0e-10;
               int max_iter=15;
               int iter=0;
               //std::cout << "Iteration #    Solution Values      Residual" << std::endl;
               //std::cout << "-----------    ---------------      --------" << std::endl;
               //std::cout << std::setw(8) << iter;
               //std::cout << std::scientific << std::setw(15) << u;
               //std::cout << std::scientific << std::setw(15) << v;
               //std::cout << std::scientific << std::setw(15) << resid << std::endl;
               while(resid>tolerance && iter<max_iter && v >= 0.0){

                  // Jacobian
                  double j11=dot3d(rho_1,rho_1)-dot3d(rho_11,dr);
                  double j12=dot3d(rho_1,rho_2)-dot3d(rho_12,dr);
                  double j21=dot3d(rho_2,rho_1)-dot3d(rho_21,dr);
                  double j22=dot3d(rho_2,rho_2)-dot3d(rho_22,dr);
                  double det=j11*j22-j12*j21;

                  // Compute increment and update newton iteration
                  // u+du = inv(Jacobian) * dr
                  u+=(j22*r1-j12*r2)/det;
                  v+=(-j21*r1+j11*r2)/det;

                  // Update residual in preparation for next iteration
                  compute_rho(u,v,rho);
                  dr[0]=xyz[0]-rho[0];
                  dr[1]=xyz[1]-rho[1];
                  dr[2]=xyz[2]-rho[2];
                  compute_rho_1(u,v,rho_1);
                  compute_rho_2(u,v,rho_2);

                  // Compute new residual
                  r1=dot3d(rho_1,dr);
                  r2=dot3d(rho_2,dr);
                  resid=sqrt(r1*r1+r2*r2);

                  // Prepare for new Jacobian
                  compute_rho_1(u,v,rho_1);
                  compute_rho_2(u,v,rho_2);
                  compute_rho_11(u,v,rho_11);
                  compute_rho_12(u,v,rho_12);
                  compute_rho_21(u,v,rho_21);
                  compute_rho_22(u,v,rho_22);

                  // print residual for this iteration
                  iter+=1;
                  //std::cout << std::setw(8) << iter;
                  //std::cout << std::scientific << std::setw(15) << u;
                  //std::cout << std::scientific << std::setw(15) << v;
                  //std::cout << std::scientific << std::setw(15) << resid << std::endl;

               }
            }
            
            // Store computed solution
            *U=u;
            *(U+1)=v;
         }

         void compute_p(double u, double *p) const {
            compute_curve(u,p);
            p[2]=-_h;
         }
         void compute_dp(double u, double *dp) const {
            compute_dcurve(u,dp);
            dp[2]=0;
         }
         void compute_ddp(double u, double *ddp) const {
            compute_dcurve(u,ddp);
            ddp[2]=0;
         }
         void compute_B(double v, double *B) const {
            B[0]=1+v*(-1+_alpha)*(2-v+2*(-1+v)*_beta);
            B[1]=B[0];
            B[2]=v;
         }
         void compute_dB(double v, double *dB) const {
            dB[0]=2*(-1+_alpha)*(1-_beta+v*(-1+2*_beta));
            dB[1]=dB[0];
            dB[2]=1;
         }
         void compute_ddB(double v, double *ddB) const {
            ddB[0]=2*(-1+_alpha)*(-1+2*_beta);
            ddB[1]=ddB[0];
            ddB[2]=0;
         }

         void compute_rho(double u, double v, double *rho) const {
            double B[]={0,0,0};
            double p[]={0,0,0};
            compute_p(u,p);
            compute_B(v,B);
            rho[0]=B[0]*p[0];
            rho[1]=B[1]*p[1];
            rho[2]=B[2]*p[2];
         }

         void compute_rho_1(double u, double v, double *rho_1) const {
            double B[]={0,0,0};
            double dp[]={0,0,0};
            compute_dp(u,dp);
            compute_B(v,B);
            rho_1[0]=B[0]*dp[0];
            rho_1[1]=B[1]*dp[1];
            rho_1[2]=B[2]*dp[2];
         }

         void compute_rho_11(double u, double v, double *rho_11) const {
            double B[]={0,0,0};
            double ddp[]={0,0,0};
            compute_ddp(u,ddp);
            compute_B(v,B);
            rho_11[0]=B[0]*ddp[0];
            rho_11[1]=B[1]*ddp[1];
            rho_11[2]=B[2]*ddp[2];
         }

         void compute_rho_12(double u, double v, double *rho_12) const {
            double dB[]={0,0,0};
            double dp[]={0,0,0};
            compute_dp(u,dp);
            compute_dB(v,dB);
            rho_12[0]=dB[0]*dp[0];
            rho_12[1]=dB[1]*dp[1];
            rho_12[2]=dB[2]*dp[2];
         }

         void compute_rho_2(double u, double v, double *rho_2) const {
            double dB[]={0,0,0};
            double p[]={0,0,0};
            compute_p(u,p);
            compute_dB(v,dB);
            rho_2[0]=dB[0]*p[0];
            rho_2[1]=dB[1]*p[1];
            rho_2[2]=dB[2]*p[2];
         }

         void compute_rho_21(double u, double v, double *rho_21) const {
            double dB[]={0,0,0};
            double dp[]={0,0,0};
            compute_dp(u,dp);
            compute_dB(v,dB);
            rho_21[0]=dB[0]*dp[0];
            rho_21[1]=dB[1]*dp[1];
            rho_21[2]=dB[2]*dp[2];
         }

         void compute_rho_22(double u, double v, double *rho_22) const {
            double ddB[]={0,0,0};
            double p[]={0,0,0};
            compute_p(u,p);
            compute_ddB(v,ddB);
            rho_22[0]=ddB[0]*p[0];
            rho_22[1]=ddB[1]*p[1];
            rho_22[2]=ddB[2]*p[2];
         }
   };

}

}

#endif
