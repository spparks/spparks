#include <cmath>
#include <stdio.h>
#include <array>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <string>
#include <sstream>
#include "am_bezier.h"
#include "fourth_order_bezier.h"

using std::pow;
using std::sqrt;
using std::vector;
using std::string;
using std::runtime_error;
using std::stringstream;

namespace bezier_pool_geometry {

   class negative_jacobian : public runtime_error {
      public:
         negative_jacobian(string& s, const array<double,3>& xyz,int iter, const array<double,2>& _u) 
            : runtime_error(s), err(s), p{xyz[0],xyz[1],xyz[2]}, u{_u[0],_u[1]}, iteration(iter) { }

            string mywhat() const {
               stringstream ss;
               ss << "\tClosest point solver encountered negative jacobian\n";
               ss << "\t\titeration = " << std::to_string(iteration) << "\n";
               ss << "\t\txyz point = " << std::scientific << std::setw(15) << p[0];
               ss << std::scientific << std::setw(15) << p[1];
               ss << std::scientific << std::setw(15) << p[2] << "\n";
               ss << "\t\tu =" << std::scientific << std::setw(15) << u[0];
               ss << std::scientific << std::setw(15) << u[1] << "\n";
               return ss.str();
            }
            array<double, 3> search_point() const { return p; }
            array<double, 2> parametric_coordinate() const { return u; }

      private:
         string err;
         array<double,3> p;
         array<double,2> u;
         int iteration;
   };

   array<double,3> compute_rho(const array<double,5>& cx, const array<double,5>& cy, const array<double,5>& cz, 
                     const array<double,4>& ry, const array<double,4>& rz, double u, double v)
   {
      // U[0]=1
      // U[1]=u
      // U[2]=pow(u,2)
      // U[3]=pow(u,3)
      // U[4]=pow(u,4)
      // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
      // doUble valUe = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
      // dc =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
      // ddc =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
      const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
      const array<double,4> V{1,pow(v,1),pow(v,2),pow(v,3)};

      double tx=0, ty=0, sz=0;
      for (int i=0;i<5;i++){
         tx+=cx[i]*U[i];
         ty+=cy[i]*U[i];
         sz+=cz[i]*U[i];
      }

      double Ty=0, Sz=0;
      for(int i=0;i<4;i++){
         Ty+=ry[i]*V[i];
         Sz+=rz[i]*V[i];
      }
      // Compute surface rho 
      array<double,3> rho{tx,ty*Ty,sz*Sz};

      return rho;

   }

   void compute_rho1(const array<double,5>& cx, const array<double,5>& cy, const array<double,5>& cz, 
                     const array<double,4>& ry, const array<double,4>& rz,  double u, double v, 
                     array<double,3>& rho, array<double,3>& rho1)
   {
      // U[0]=1
      // U[1]=u
      // U[2]=pow(u,2)
      // U[3]=pow(u,3)
      // U[4]=pow(u,4)
      // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
      // doUble valUe = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
      // dc =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
      // ddc =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
      const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
      const array<double,4> V{1,pow(v,1),pow(v,2),pow(v,3)};

      double tx=0, ty=0, sz=0;
      for (int i=0;i<5;i++){
         tx+=cx[i]*U[i];
         ty+=cy[i]*U[i];
         sz+=cz[i]*U[i];
      }
      // 1st derivatives
      double dtx=0, dty=0, dsz=0;
      for (int i=1;i<5;i++){
         dtx+=i*cx[i]*U[i-1];
         dty+=i*cy[i]*U[i-1];
         dsz+=i*cz[i]*U[i-1];
      }

      double Ty=0, Sz=0;
      for(int i=0;i<4;i++){
         Ty+=ry[i]*V[i];
         Sz+=rz[i]*V[i];
      }

      double dTy=0, dSz=0;
      for(int i=0;i<4;i++){
         dTy+=ry[i]*V[i];
         dSz+=rz[i]*V[i];
      }
      // Compute surface rho
      // Compute surface rho1; derivative wrt u
      rho[0]=tx;
      rho[1]=ty*Ty;
      rho[2]=sz*Sz;
      rho1[0]=dtx;
      rho1[1]=dty*Ty;
      rho1[2]=dsz*Sz;

   }

   void compute_rho2(const array<double,5>& cx, const array<double,5>& cy, const array<double,5>& cz, 
                     const array<double,4>& ry, const array<double,4>& rz, double u, double v, 
                     array<double,3>& rho, array<double,3>& rho2)
   {
      // U[0]=1
      // U[1]=u
      // U[2]=pow(u,2)
      // U[3]=pow(u,3)
      // U[4]=pow(u,4)
      // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
      // doUble valUe = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
      // dc =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
      // ddc =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
      const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
      const array<double,4> V{1,pow(v,1),pow(v,2),pow(v,3)};

      double tx=0, ty=0, sz=0;
      for (int i=0;i<5;i++){
         tx+=cx[i]*U[i];
         ty+=cy[i]*U[i];
         sz+=cz[i]*U[i];
      }

      double Ty=0, Sz=0;
      for(int i=0;i<4;i++){
         Ty+=ry[i]*V[i];
         Sz+=rz[i]*V[i];
      }

      double dTy=0, dSz=0;
      for(int i=1;i<4;i++){
         dTy+=i*ry[i]*V[i-1];
         dSz+=i*rz[i]*V[i-1];
      }

      // Compute surface rho 
      // Compute surface rho2; derivative wrt v
      rho[0]=tx;
      rho[1]=ty*Ty;
      rho[2]=sz*Sz;
      rho2[0]=0;
      rho2[1]=ty*dTy;
      rho2[2]=sz*dSz;
   }

namespace closest_point_projection {

   inline double dot(const array<double,2>& A, const array<double,2>& B){

      return (A[0]*B[0]+A[1]*B[1]);

   }

   inline double dot(const array<double,3>& A, const array<double,3>& B){

      return (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]);

   }

   void compute_linear_system(
         const array<double,5>& cx,
         const array<double,5>& cy,
         const array<double,5>& cz,
         const array<double,4>& ry, const array<double,4>& rz,
         const array<double,3>& XYZ,
         double u, double v, 
         array<array<double,2>,2>& J, array<double,2>& rhs)
   {
      // U[0]=1
      // U[1]=u
      // U[2]=pow(u,2)
      // U[3]=pow(u,3)
      // U[4]=pow(u,4)
      // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
      // doUble valUe = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
      // dc =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
      // ddc =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
      const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
      const array<double,4> V{1,pow(v,1),pow(v,2),pow(v,3)};

      double tx=0, ty=0, sz=0;
      for (int i=0;i<5;i++){
         tx+=cx[i]*U[i];
         ty+=cy[i]*U[i];
         sz+=cz[i]*U[i];
      }
      // 1st derivatives
      double dtx=0, dty=0, dsz=0;
      for (int i=1;i<5;i++){
         dtx+=i*cx[i]*U[i-1];
         dty+=i*cy[i]*U[i-1];
         dsz+=i*cz[i]*U[i-1];
      }
      // 2nd derivatives
      double ddtx=0, ddty=0, ddsz=0;
      for (int i=2;i<5;i++){
         ddtx+=i*(i-1)*cx[i]*U[i-2];
         ddty+=i*(i-1)*cy[i]*U[i-2];
         ddsz+=i*(i-1)*cz[i]*U[i-2];
      }


      double Ty=0, Sz=0;
      for(int i=0;i<4;i++){
         Ty+=ry[i]*V[i];
         Sz+=rz[i]*V[i];
      }

      // 1st derivatives
      double dTy=0, dSz=0;
      for(int i=1;i<4;i++){
         dTy+=i*ry[i]*V[i-1];
         dSz+=i*rz[i]*V[i-1];
      }

      // 2nd derivatives
      double ddTy=0, ddSz=0;
      for(int i=2;i<4;i++){
         ddTy+=i*(i-1)*ry[i]*V[i-2];
         ddSz+=i*(i-1)*rz[i]*V[i-2];
      }

      // Compute surface rho 
      const array<double,3> rho{tx,ty*Ty,sz*Sz};

      // Compute surface rho1; derivative wrt u
      const array<double,3> rho1{dtx,dty*Ty,dsz*Sz};

      // Compute surface rho11; 2 derivatives wrt u
      const array<double,3> rho11{ddtx,ddty*Ty,ddsz*Sz};

      // Following two derivatives are related and similar
      // Compute surface rho2; derivative wrt v
      const array<double,3> rho2{0,ty*dTy,sz*dSz};

      // Compute surface rho12; derivative wrt u derivative wrt v
      const array<double,3> rho12{0,dty*dTy,dsz*dSz};

      // Compute surface rho22; 2 derivatives wrt v
      const array<double,3> rho22{0,ty*ddTy,sz*ddSz};

      //printf("AM_Bezier::print_surface() u=%10.4f, v=%10.4f\n",u,v);
      //printf("\trhox=%10.4f, rhoy=%10.4f, rhoz=%10.4f\n",rho[0],rho[1],rho[2]);
      //printf("\trho1x=%10.4f, rho1y=%10.4f, rho1z=%10.4f\n",rho1[0],rho1[1],rho1[2]);
      //printf("\trho11x=%10.4f, rho11y=%10.4f, rho11z=%10.4f\n",rho11[0],rho11[1],rho11[2]);
      //printf("\trho12x=%10.4f, rho12y=%10.4f, rho12z=%10.4f\n",rho12[0],rho12[1],rho12[2]);
      //printf("\trho2x=%10.4f, rho2y=%10.4f, rho2z=%10.4f\n",rho2[0],rho2[1],rho2[2]);
      //printf("\trho22x=%10.4f, rho22y=%10.4f, rho22z=%10.4f\n",rho22[0],rho22[1],rho22[2]);

      // Compute 2x2 matrix of linear system
      // right hand side vector rhs
      const array<double,3> dr{XYZ[0]-rho[0],XYZ[1]-rho[1],XYZ[2]-rho[2]};
      rhs[0]=-dot(rho1,dr);
      rhs[1]=-dot(rho2,dr);

      // 2x2 matrix
      J[0][0]=dot(rho11,dr)-dot(rho1,rho1);
      J[0][1]=dot(rho12,dr)-dot(rho1,rho2);
      J[1][0]=J[0][1];
      J[1][1]=dot(rho22,dr)-dot(rho2,rho2);
   }

   array<double,2>
   predict(const array<double,5>& cx,  const array<double,5>& cy, 
             const array<double,5>& cz,  const array<double,4>& ry, const array<double,4>& rz, 
             const array<double,3>& xyz, double umax_width, double umax_depth)
   {

      using fourth_order_bezier_polynomial::compute;
      using fourth_order_bezier_polynomial::cpp_predictor_2d;
      using fourth_order_bezier_polynomial::solve;
      using bezier_pool_geometry::compute_rho2;
      array<double,2> u0;
      const double tolerance=1.0e-15;
      double xm;
      {
         // this computes quick estimate of u
         // Identify quadrant of input spatial point (x,y)
         {
            // Compute x-coordinate at maximum pool width
            xm=compute(cx,umax_width);
         }
         // Starting value for bounds on parametric coordinate
         double umin, umax;
         double x=xyz[0];
         double y=xyz[1];
         double z=xyz[2];
         if (x>xm){
            umin=umax_width;
            umax=1.0;
         }
         else {
            umin=0;
            umax=umax_width;
         }
         // u estimate
         u0[0]=0.5*(umin+umax);

         // override for points on spine at top surface
         const double head=compute(cx,1.0);
         const double tail=compute(cx,0.0);
         if(y<tolerance){
            array<double,2> xz{xyz[0],xyz[2]};
            double _u=cpp_predictor_2d(cx,cz,xz,umax_depth,.001);
            u0[0]=_u;
            if(fabs(z)<tolerance){
               // points at top surface
               u0[1]=0;
               return u0;
            } else{
               // points on spine but below surface
               u0[1]=1.0;
               return u0;
            }
         }

         {
            // if x between head and tail,
            // then can solve u exactly on top curve
            if(x<head && x>tail){
               const int max_iter=15;
               u0[0]=solve(cx,u0[0],x,tolerance,max_iter);
            }
         }

         // at top surface v=0
         if(abs(z)<tolerance){
            u0[1]=0.0;
         } else {
            // Bound parametic coordinate and compute (umin, umax) of bounds
            double vmin=0,vmax=0.5;
            double dv=.1;
            double vn=0.5*(vmin+vmax);
            while (vmax-vmin > dv){
               array<double,3> rho;
               array<double,3> rho2;
               compute_rho2(cx,cy,cz,ry,rz,u0[0],vn,rho,rho2);
               array<double,3> dr={x-rho[0],y-rho[1],z-rho[2]};
               double dt=dot(dr,rho2);
               if (dt>0){
                   // Point below un
                   vmin=vn;
               }
               else{
                   // Point above un
                   vmax=vn;
               }
               vn=0.5*(vmin+vmax);
            }
            u0[1]=vn;
         }
      }

      return u0;

   }

   array<double, 2> compute_closest_point(
                  const array<double,5>& cx,
                  const array<double,5>& cy,
                  const array<double,5>& cz,
                  const array<double,4>& ry, 
                  const array<double,4>& rz,
                  const array<double,3>& xyz, 
                  const array<double,2>& u0) 
   {
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
       * u0: 2 component vector storing initial estimate of closest point projection.
       *    
       * Return value
       * -------
       *  Computed closest point U in parametric coordinates of surface.
       * 
       */


      // Initial solution guess/estimate
      double u=u0[0];
      double v=u0[1];
      // Jacobian 2x2 matrix
      array<array<double,2>,2> J;
      // residual vector
      array<double,2> rhs{0,0};
      double ref=sqrt(dot(xyz,xyz));
      compute_linear_system(cx,cy,cz,ry,rz,xyz,u,v,J,rhs);
      double resid=sqrt(dot(rhs,rhs));
      {

         double tolerance=1.0e-10;
         int max_iter=15;
         int iter=0;
         //std::cout << "Iteration #    Solution Values      Residual" << std::endl;
         //std::cout << "-----------    ---------------      --------" << std::endl;
         double j11, j12, j21, j22, det, r1, r2;
         while((resid)/ref>tolerance && iter<max_iter) {

            // Jacobian 
            j11=J[0][0];
            j12=J[0][1];
            j21=J[1][0];
            j22=J[1][1];
            det=j11*j22-j12*j21;

            if(det<0){
               //throw std::runtime_error("AM_Bezier.cpp; compute_closest point; 'negative jacobian'\n");
               string s("Closest point projection solver encountered negative jacobian");
               throw negative_jacobian(s,xyz,iter,array<double,2>{u,v});
            }

            // rhs
            r1=rhs[0];
            r2=rhs[1];

            // iteration print
            //std::cout << std::setw(8) << iter;
            //std::cout << std::scientific << std::setw(15) << u;
            //std::cout << std::scientific << std::setw(15) << v;
            //std::cout << std::scientific << std::setw(15) << det;
            //std::cout << std::scientific << std::setw(15) << resid/ref << std::endl;

            // Compute increment and update newton iteration
            // u+du = inv(Jacobian) * dr
            //u+=(j22*r1-j12*r2)/det;
            double du=(j22*r1-j12*r2)/det;
            //std::cout << "\t du = " << std::scientific << std::setw(15) << du;
            // HACKs to prevent overshoot 
            // Better predictor may solve this problem
            while(u+du>1.0 || u+du<0) du*=.5;
            u+=du;
            // HACKs to prevent overshoot 
            // Better predictor may solve this problem
            double dv=(-j21*r1+j11*r2)/det;
            //std::cout << "\t dv = " << std::scientific << std::setw(15) << dv << std::endl;
            //v+=(-j21*r1+j11*r2)/det;
            while(v+dv>1.0 || v+dv<0.0) dv*=.5;
            v+=dv;

            // Update residual and Jacobian in preparation for next iteration
            compute_linear_system(cx,cy,cz,ry,rz,xyz,u,v,J,rhs);
            resid=sqrt(dot(rhs,rhs));

            // print residual for this iteration
            iter+=1;

         }
         //std::cout << std::setw(8) << iter;
         //std::cout << std::scientific << std::setw(15) << u;
         //std::cout << std::scientific << std::setw(15) << v;
         //std::cout << std::scientific << std::setw(15) << det;
         //std::cout << std::scientific << std::setw(15) << resid/ref << std::endl;

      }
      
      // Solution is return value;
      array<double,2> U;
      U[0]=u;
      U[1]=v;
      return U;
   }

} // closest_point_projection

} // bezier_pool_geometry

AM_Bezier::
AM_Bezier() :
	xt{0,0,0,0,0}, 
	yt{0,0,0,0,0}, 
	zs{0,0,0,0,0}, 
	ry{0,0,0,0}, 
	rz{0,0,0,0}, 
   width(1.0),umax_width(-1.0),length(1.0),
   depth(1.0), umax_depth(-1.0),
   bb_center{0,0,0},bb_dimensions{0,0,0},
   search_points(0),search_coordinates(0),
   tree(my_tree::get_tree(nullptr,0))
{ }

AM_Bezier::
~AM_Bezier(){ }

double AM_Bezier::
compute_pool_length() const {

   using fourth_order_bezier_polynomial::compute;
   // computes current pool length
   // based upon precomputed cx coefficients
   //@u=0.0
   double x0=compute(cx,0.0);
   //@u=1.0
   double x1=compute(cx,1.0);
   // computing private variable
   double length=std::fabs(x1-x0);

   return length;
}

AM_Bezier::
AM_Bezier(
          const array<double,5>& x,
          const array<double,5>& y,
          const array<double,5>& z,
          double betay,
          double betaz,
          double width,
          double depth
         ) :
	xt{x[0],x[1],x[2],x[3],x[4]}, 
	yt{y[0],y[1],y[2],y[3],y[4]}, 
	zs{z[0],z[1],z[2],z[3],z[4]}, 
	cx{0.0,0.0,0.0,0.0,0.0}, 
	cy{0.0,0.0,0.0,0.0,0.0}, 
	cz{0.0,0.0,0.0,0.0,0.0}, 
	ry{1,0,(-3+3*betay),(2-3*betay)},
	rz{0,3*betaz,(3-6*betaz),(-2+3*betaz)},
   width(width),length(1.0),depth(depth),
   umax_width(-1.0), umax_depth(-1.0),
   bb_center{0,0,0},bb_dimensions{0,0,0},
   search_points(0),search_coordinates(0),
   tree(my_tree::get_tree(nullptr,0))
{

   using fourth_order_bezier_polynomial::compute_coefficients;
   using fourth_order_bezier_polynomial::scale_coordinates_for_curve_height;

   // updates coordinates xt, yt coordinates for HALF width
   std::tie(umax_width,xt,yt)=scale_coordinates_for_curve_height(x,y,0.5*width);

   // updates zt for depth
   // be careful here; do not re-define x-coordinates based on depth scaling;
   // only use scaled z-coordinates
   array<double,5> unused;
   std::tie(umax_depth,unused,zs)=scale_coordinates_for_curve_height(x,z,depth);
   cx=compute_coefficients(xt);
   cy=compute_coefficients(yt);
   cz=compute_coefficients(zs);
   length=compute_pool_length();
   double xfront=xt[4];
   double xtail=xt[0];
   double xc=0.5*(xfront+xtail);
   double yc=0;
   double zc=-0.5*depth;
   double hx=length;
   double hy=width;
   double hz=depth;
   bb_center[0]=xc;
   bb_center[1]=yc;
   bb_center[2]=zc;
   bb_dimensions[0]=hx;
   bb_dimensions[1]=hy;
   bb_dimensions[2]=hz;

   // create kdtree points
   {
      // for kdtree search
      // create points forward of max width
      using bezier_pool_geometry::compute_rho;
      //double u0=umax_width, u1=0.99, v0=0, v1=1.0;
      double u0=.01, u1=0.99, v0=0, v1=1.0;
      int nu=100, nv=100;
      double du=(u1-u0)/nu;
      double dv=(v1-v0)/nv;
      int num_points=nu*nv;
      int count(0);
      for(double u=u0;u<=u1;u+=du){
         for(double v=v0;v<=v1;v+=dv){
            array<double,3> rho = compute_rho(cx, cy, cz, ry, rz, u, v);
            double x=rho[0];
            double y=rho[1];
            double z=rho[2];
            search_points.push_back(x);
            search_points.push_back(y);
            search_points.push_back(z);
            search_coordinates.push_back(array<double,2>{u,v});
            //std::cout << std::setw(8) << "point i: " << count;
            //std::cout << std::scientific << std::setw(15) << u;
            //std::cout << std::scientific << std::setw(15) << v;
            //std::cout << std::scientific << std::setw(15) << x;
            //std::cout << std::scientific << std::setw(15) << y;
            //std::cout << std::scientific << std::setw(15) << z << std::endl;
            count+=1;
         }
      }
      //printf("AM_Bezier creating %d search points \n",count);
      //printf("AM_Bezier constructor; num_points = %d\n",dum);
      //for(int i=0;i<num_points;i++){
      //   double x=search_points[3*i];
      //   double y=search_points[3*i+1];
      //   double z=search_points[3*i+2];
      //   std::cout << std::setw(8) << "point i: " << i;
      //   std::cout << std::scientific << std::setw(15) << x;
      //   std::cout << std::scientific << std::setw(15) << y;
      //   std::cout << std::scientific << std::setw(15) << z << std::endl;
      //}

      tree=my_tree::get_tree(search_points.data(),num_points);
   }

   /*
   for(int i=0;i<search_points.size()/3;i++){
      vector<double> p{search_points[3*i],search_points[3*i+1],search_points[3*i+2]};
      int id=tree->nearest_neighbor_search(p);
      std::cout << std::setw(8) << "search point i: " << i << "; found id = " << id;
      std::cout << std::scientific << std::setw(15) << p[0];
      std::cout << std::scientific << std::setw(15) << p[1];
      std::cout << std::scientific << std::setw(15) << p[2] << std::endl;
   }
   */

}


void 
AM_Bezier::
print_rho(double u, double v) const 
{


}

bool
AM_Bezier::
is_outside_pool(const array<double,3>& xyz) const {

   // Until more accurate method is developed, 
   //    use pool bounding box
   double xmin=bb_center[0]-bb_dimensions[0]/2;
   double xmax=bb_center[0]+bb_dimensions[0]/2;
   double ymin=bb_center[1]-bb_dimensions[1]/2;
   double ymax=bb_center[1]+bb_dimensions[1]/2;
   double zmin=bb_center[2]-bb_dimensions[2]/2;
   double zmax=bb_center[2]+bb_dimensions[2]/2;
   double x=xyz[0];
   double y=xyz[1];
   double z=xyz[2];

   bool outside=false;
   if(x<xmin || x>xmax)
      outside=true;
   else if(y<ymin || y>ymax)
      outside=true;
   else if(z<zmin || z>zmax)
      outside=true;

   //printf("AM_Bezier::is_outside_pool = %d\n",outside);
   return outside;
}

bool
AM_Bezier::
inside_outside(const array<double,3>& xyz, const double tolerance) const {
   using fourth_order_bezier_polynomial::cpp_predictor_2d;
   using fourth_order_bezier_polynomial::solve;
   using bezier_pool_geometry::compute_rho;
   using bezier_pool_geometry::compute_rho2;
   using bezier_pool_geometry::closest_point_projection::dot;

   //check if inside boundary box
   //returns true if outside and false if inside
   bool outside=is_outside_pool(xyz);
   if (!outside){ 
      double _u=0.5;
      double _v=0.5;
      {
         array<double,2> xy{xyz[0],xyz[1]};
         double u0=cpp_predictor_2d(cx,cy,xy,umax_width,.001);
         const int max_iter=15;
         _u=solve(cx,u0,xy[0],tolerance,max_iter);
      }

      //get z limits for current cross-section (for our given _u)
      array<double,3> rho = compute_rho(cx, cy, cz, ry, rz, _u, 1);
      double zmin = rho[2];
      double ymin = 0;
      rho = compute_rho(cx, cy, cz, ry, rz, _u, 0);
      double zmax = rho[2];
      double ymax=rho[1];

      if (xyz[2] <= zmax && xyz[2] >= zmin && xyz[1] <= ymax && xyz[1] >= ymin){
         //solve for the v value corresponding to xyz[2]   
         const int max_iter=25;
         int iter=0;
         array<double,3> rho;
         array<double,3> rho2;
         compute_rho2(cx, cy, cz, ry, rz, _u, _v, rho, rho2);
         double f=-(rho[2]-xyz[2]);
         double j=rho2[2];

         double resid=std::fabs(f);
         //std::cout << "Iteration #    v               dy          dz                 f          Residual" << std::endl;
         //std::cout << "-----------  ---------       -------      ---------       ---------     ---------" << std::endl;
         while(resid>tolerance && ++iter<=max_iter){
            //std::cout << std::setw(8) << iter;
            //std::cout << std::scientific << std::setw(15) << _v;
            //std::cout << std::scientific << std::setw(15) << rho2[1];
            //std::cout << std::scientific << std::setw(15) << j;
            //std::cout << std::scientific << std::setw(15) << f;
            //std::cout << std::scientific << std::setw(15) << resid << std::endl;
            double dv = f/j;
            while(_v+dv >= 1 || _v+dv <= 0) dv*=0.5;
            _v+=dv;
            compute_rho2(cx, cy, cz, ry, rz, _u, _v, rho, rho2);
            f=-(rho[2]-xyz[2]);
            j=rho2[2];
            resid=std::fabs(f);
         }

         rho = compute_rho(cx, cy, cz, ry, rz, _u, _v);

         if (xyz[1]<rho[1]){
            return 0; // return 0 for inside
         }
      }
   }
   return 1; // return 1 for outside
}

double 
AM_Bezier::
distance(const array<double,3>& _xyz) const {

   // xyz is position of point in space relative to 
   //   reference point {0,0,0} of pool
   // Depending of x-coordinates of control points, 
   //   reference point should be along centerline 
   //   of pool somewhere between front and tail; 
   //   Probably best for it to be somewhere towards 
   //   the middle
   // Y-coordinate of reference point is on centerline 
   // Reference z is at top which is taken as zero 
   //   in the local pool coordinate system.
   // Pool travels along centerline which x-axis

   // Search point which is straight away from maximum width
   using bezier_pool_geometry::closest_point_projection::predict;
   using bezier_pool_geometry::closest_point_projection::compute_closest_point;
   using bezier_pool_geometry::compute_rho;
   using bezier_pool_geometry::closest_point_projection::dot;

   // Need to flip point to +y
   // Need to determine if point is inside or out
   // Identify points at surface z=0 or |z|<z_top
   //   * points at surface, do cpp on top curve
   //   * calculate distance 
   //   * return distance
   // If point not on top surface continue
   // If point on pool centerline |y|<tolerance
   //   * point is on spine
   //   * do cpp on spine curve
   //   * calculate distance
   //   * return distance
   // If point is not on surface or spine
   //   * need predictor
   //   * do cpp on surface

   // flip y-coordinate to positive 
   array<double,3> xyz{_xyz[0],fabs(_xyz[1]),_xyz[2]};
   double d=-1.0;
   {
      using fourth_order_bezier_polynomial::closest_point;
      using fourth_order_bezier_polynomial::dot;
      const double tolerance=1.0e-15;
      // check if top plane
      if(abs(xyz[2])<tolerance){
         array<double,2> xy{xyz[0],xyz[1]};
         const int max_iter=15;
         double unused;
         std::tie(unused,d)=closest_point(cx,cy,xy,umax_width,tolerance,max_iter);
         return d;
      }
      // check if spine
      // BIG NOTE
      // Notice dot product function is passed to cpp below
      // This is because closest point project depends up 
      //    CCW or CW path of curve
      // cpp in xz plane uses CCW path which is right handed
      // cpp in xz plane uses CW path which is LEFT handed (default case)
      if(xyz[1]<tolerance){
         array<double,2> xz{xyz[0],xyz[2]};
         const int max_iter=15;
         double unused;
         std::tie(unused,d)=closest_point(cx,cz,xz,umax_depth,tolerance,max_iter,dot);
         return d;
      }
   }

   // test point for inside or outside of pool
   //bool outside=is_outside_pool(xyz);
   //bool outside=is_outside_pool(xyz);
   bool outside=inside_outside(xyz,1.0E-10);

   // closest point projection and distance to pool surface
   if(outside){
      array<double,2> u;
      array<double,2> u0=predict(cx,cy,cz,ry,rz,xyz,umax_width,umax_depth);
      try {
         u=compute_closest_point(cx,cy,cz,ry,rz,xyz,u0);
      } catch(bezier_pool_geometry::negative_jacobian& e){
         //std::cout << e.mywhat() << std::endl;
         int u0id=tree.nearest_neighbor_search(vector<double>{xyz[0],xyz[1],xyz[2]});
         u0=search_coordinates[u0id];
         //printf("\t\tTrying kdtree search value u0={%10.6f, %10.6f}\n",u0[0],u0[1]);
         try{
            u=compute_closest_point(cx,cy,cz,ry,rz,xyz,u0);
         } catch(...){
            // Give up on cpp; use the kdtree result
            u=u0;
         }
      } 
      array<double,3> rho=compute_rho(cx,cy,cz,ry,rz,u[0],u[1]);
      array<double,3> dr{xyz[0]-rho[0],xyz[1]-rho[1],xyz[2]-rho[2]};
      //printf("\tpredicted solution u0={%10.6f, %10.6f}\n",u0[0],u0[1]);
      //printf("\tsolution u={%10.6f, %10.6f}\n",u[0],u[1]);
      //printf("\tsolution rho={%10.6f, %10.6f, %10.6f}\n",rho[0],rho[1],rho[2]);
      //printf("\tsolution dr={%10.6f, %10.6f, %10.6f}\n",dr[0],dr[1],dr[2]);
      d=sqrt(dot(dr,dr));
   }
   return d;
}

void AM_Bezier::print() const
{
   double l = get_pool_length();
   double d = get_pool_depth();
   double w = get_pool_width();
   printf("AM_Bezier::print() pool width=%16.13f @ umax_width=%16.13f\n",w,umax_width);
   printf("AM_Bezier::print() pool depth=%16.13f @ umax_depth=%16.13f\n",d,umax_depth);
   printf("AM_Bezier::print() pool length=%16.13f\n",l);
   printf("AM_Bezier::print() length/width=%16.13f\n",l/w);
   printf("AM_Bezier::print() length/depth=%16.13f\n",l/d);
   double bbx=bb_center[0];
   double bby=bb_center[1];
   double bbz=bb_center[2];
   double xmin=bb_center[0]-bb_dimensions[0]/2;
   double ymin=bb_center[1]-bb_dimensions[1]/2;
   double zmin=bb_center[2]-bb_dimensions[2]/2;
   double xmax=bb_center[0]+bb_dimensions[0]/2;
   double ymax=bb_center[1]+bb_dimensions[1]/2;
   double zmax=bb_center[2]+bb_dimensions[2]/2;
   printf("AM_Bezier::print() bounding box center x,y,z=[%8.2f,%8.2f,%8.2f]\n",bbx,bby,bbz);
   printf("AM_Bezier::print() bounding box x=[%8.2f,%8.2f], y=[%8.2f,%8.2f],z=[%8.2f,%8.2f]\n",xmin,xmax,ymin,ymax,zmin,zmax);
   //const double *xo=this->get_pool_reference_center();
   //printf("AM_Bezier::print() pool reference center=%10.4f, %10.4f, %10.4f\n",xo[0],xo[1],xo[2]);
}
