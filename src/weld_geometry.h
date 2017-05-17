#ifndef SPK_WELD_GEOMETRY_H
#define SPK_WELD_GEOMETRY_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include "pool_shape.h"

namespace weld {

namespace pool_shape {

   class ParameterizedEllipse {

      public:
         ParameterizedEllipse(double _x_width, double _y_width): 
            a(_x_width), b(_y_width), a2(a*a), b2(b*b) { } 

         void get_coordinates(double theta, double *xy) const {
           xy[0]=a*cos(theta); 
           xy[1]=b*sin(theta); 
         }

         double closest_point(double x, double y) const {
            double theta;

            {
               /*
                * Solve for 't' that minimizes distance
                */
               int max_iter=20;
               double tolerance=1.0e-15;
               int iter=0;
               double t=0;
               double r=f(x,y,t);
               while (r>tolerance && iter<max_iter){
                  t=t-r/df(x,y,t);
                  r=f(x,y,t);
                  iter+=1;
               }
               // spatial coordinates of point on ellipse which minimizes distance
               double x0=a2*x/(a2+t);
               double y0=b2*y/(b2+t);

               // Invert for parametric coordinate 'theta'
               double ac;
               if(x0/a+1.0<0.0){ ac=M_PI; } 
               else if(x0/a-1.0>0.0){ ac=0; } 
               else {ac=acos(x0/a);}

               if (y0>=0.0){ 
                  // Handles quadrants I and II
                  theta=ac;
               }
               else {
                  // Handles quadrants III and IV
                  theta=2*M_PI-ac;
               }

            }

            return theta;
         }


      private:
        double f(double x, double y, double t) const {
           double ax=a*x;
           double by=b*y;
           double c1=ax/(a2+t);
           double c2=by/(b2+t);
           return pow(c1,2)+pow(c2,2)-1.0;
        }

        double df(double x, double y, double t) const {
           double ax=a*x;
           double by=b*y;
           double c1=-2.0*pow(ax,2)/pow(a2+t,3);
           double c2=-2.0*pow(by,2)/pow(b2+t,3);
           return c1+c2;
        }


      private:
         const double a, b;
         const double a2, b2;

   };

   class EllipticBezier : public PoolShape {
      private:
         const double a, b, T, alpha, beta;

      public:
         EllipticBezier 
            (
              double pool_width, double pool_length, double plate_thickness,
              double scale_param, double interpolate_param
            ):
            a(pool_width/2), b(pool_length/2.0), T(plate_thickness),
            alpha(scale_param), beta(interpolate_param)
         { } 

         virtual ~EllipticBezier() {}
         
         bool is_inside(const double *xyz) const {
            double x=*xyz;
            double y=*(xyz+1);
            double z=*(xyz+2);
            double x2=x*x;
            double y2=y*y;
            double z2=z*z;
            double T2=T*T;
            double c=(T2+2*T*z*(-1+alpha)*(-1+beta)+z2*(-1+alpha)*(-1+2*beta))/T2;
            double _a=a*c;
            double _b=b*c;
            double _a2=_a*_a;
            double _b2=_b*_b;
            return x2/_a2 + y2/_b2 < 1.0;
         }

         virtual double distance(const double *xyz) const {
            /**
             * For points 'xyz' outside of pool, 
             * computes closest point projection to 
             * 'EllipticBezier' surface; Using this 
             * projection the distance to the surface 
             * is returned.
             */

            // Must be outside of pool
            double d=-1.0;
            if(is_inside(xyz)) return d;

            // Parameteric coordinates for closest point project
            double u0, v0;
            {
               // Compute initial guess (u0,v0) at parametric coordinates
               double x=*xyz;
               double y=*(xyz+1);
               double z=*(xyz+2);
               double z2=z*z;
               double T2=T*T;
               double c=(T2+2*T*z*(-1+alpha)*(-1+beta)+z2*(-1+alpha)*(-1+2*beta))/T2;
               double _a=a*c;
               double _b=b*c;
               ParameterizedEllipse pe(_a,_b);
               u0=pe.closest_point(x,y);
               v0=-z/T;
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
            d=sqrt(dot(dr,dr));
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
               double r1=dot(rho_1,dr);
               double r2=dot(rho_2,dr);
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
                  double j11=dot(rho_1,rho_1)-dot(rho_11,dr);
                  double j12=dot(rho_1,rho_2)-dot(rho_12,dr);
                  double j21=dot(rho_2,rho_1)-dot(rho_21,dr);
                  double j22=dot(rho_2,rho_2)-dot(rho_22,dr);
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
                  r1=dot(rho_1,dr);
                  r2=dot(rho_2,dr);
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
            p[0]=a * cos(u);
            p[1]=b * sin(u);
            p[2]=-T;
         }
         void compute_dp(double u, double *dp) const {
            dp[0]=-a * sin(u);
            dp[1]= b * cos(u);
            dp[2]=0;
         }
         void compute_ddp(double u, double *ddp) const {
            ddp[0]=-a * cos(u);
            ddp[1]=-b * sin(u);
            ddp[2]=0;
         }
         void compute_B(double v, double *B) const {
            B[0]=1+v*(-1+alpha)*(2-v+2*(-1+v)*beta);
            B[1]=B[0];
            B[2]=v;
         }
         void compute_dB(double v, double *dB) const {
            dB[0]=2*(-1+alpha)*(1-beta+v*(-1+2*beta));
            dB[1]=dB[0];
            dB[2]=1;
         }
         void compute_ddB(double v, double *ddB) const {
            ddB[0]=2*(-1+alpha)*(-1+2*beta);
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

         double dot(const double* y1, const double *y2) const { 
            return y1[0]*y2[0]+y1[1]*y2[1]+y1[2]*y2[2];
         }

   };


}

}

#endif
