#include "fourth_order_bezier.h"
#include <cmath>
#include <cstdio>

using std::pow;
using std::sqrt;
using std::make_tuple;


namespace fourth_order_bezier_polynomial {


double cw_path_dot(const array<double,2>& A, const array<double,2>& B){

   return -(A[0]*B[0]+A[1]*B[1]);

}

double dot(const array<double,2>& A, const array<double,2>& B){

   return (A[0]*B[0]+A[1]*B[1]);

}

array<double,5> compute_coefficients(const array<double,5>& x) {
   // x[5] control points for 4th order bezier polynomial
   // double value = c0 + c1 * u + c2 * pow(u,2) + c3 * pow(u,3) + c4 * pow(u,4);
   double x0=x[0];
   double x1=x[1];
   double x2=x[2];
   double x3=x[3];
   double x4=x[4];

   array<double,5> c;
   c[0]=x0;
   c[1]=(-4*x0 + 4*x1);
   c[2]=(6*x0 - 12*x1 + 6*x2) ;
   c[3]=(-4*x0 + 12*x1 - 12*x2 + 4*x3);
   c[4]=(x0 - 4*x1 + 6*x2 - 4*x3 + x4);

   return c;
}

double compute(const array<double,5>& coefficients, double u) {
   // coefficients[5] derived from precomputed control points
   // double value = c0 + c1 * u + c2 * pow(u,2) + c3 * pow(u,3) + c4 * pow(u,4);
   double u2 = pow(u,2);
   double u3 = pow(u,3);
   double u4 = pow(u,4);
   double c0=coefficients[0];
   double c1=coefficients[1];
   double c2=coefficients[2];
   double c3=coefficients[3];
   double c4=coefficients[4];
   double value = c0 + c1 * u + c2 * u2 + c3 * u3 + c4 * u4;
   return value;
}

void compute(const array<double,5>& coefficients, double u, array<double,3>& c) {
   // coefficients[5] derived from precomputed control points
   // double value = c0 + c1 * u + c2 * pow(u,2) + c3 * pow(u,3) + c4 * pow(u,4);
   double u2 = pow(u,2);
   double u3 = pow(u,3);
   double u4 = pow(u,4);
   double c0=coefficients[0];
   double c1=coefficients[1];
   double c2=coefficients[2];
   double c3=coefficients[3];
   double c4=coefficients[4];
   // curve
   c[0] = c0 + c1 * u + c2 * u2 +   c3 * u3 +    c4 * u4;
   // 1st derivative of curve
   c[1] =      c1   + 2*c2 * u  + 3*c3 * u2 +  4*c4 * u3;
   // 2nd derivative of curve
   c[2] =             2*c2      + 6*c3 * u  + 12*c4 * u2;
}


void curve(const array<double,5>& cx, const array<double,5>& cy, double u, array<double,2>& v) {
   // cx, cy precomputed coefficients
   // U[0]=1
   // U[1]=u
   // U[2]=pow(u,2)
   // U[3]=pow(u,3)
   // U[4]=pow(u,4)
   // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
   // double f = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
   // df =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
   // ddf =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
   const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
   // initialize 
   for(int i=0;i<2;i++){
      v[i]=0;
   }
   // compute curve
   for (int i=0;i<5;i++){
      v[0]+=cx[i]*U[i];
   }
}

void curve(const array<double,5>& cx, const array<double,5>& cy, double u, 
           array<double,2>& v, array<double,2>& dv) {
   // cx, cy precomputed coefficients
   // U[0]=1
   // U[1]=u
   // U[2]=pow(u,2)
   // U[3]=pow(u,3)
   // U[4]=pow(u,4)
   // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
   // double f = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
   // df =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
   // ddf =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
   const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
   // initialize 
   for(int i=0;i<2;i++){
      v[i]=0;
      dv[i]=0;
   }
   // compute curve
   for (int i=0;i<5;i++){
      v[0]+=cx[i]*U[i];
      v[1]+=cy[i]*U[i];
   }
   // compute 1st derivative of curve wrt u
   for (int i=1;i<5;i++){
      dv[0]+=i*cx[i]*U[i-1];
      dv[1]+=i*cy[i]*U[i-1];
   }
}

void curve(const array<double,5>& cx, const array<double,5>& cy, double u, 
           array<double,2>& v, array<double,2>& dv, array<double,2>& ddv) {
   // cx, cy precomputed coefficients
   // U[0]=1
   // U[1]=u
   // U[2]=pow(u,2)
   // U[3]=pow(u,3)
   // U[4]=pow(u,4)
   // coefficients cx[5], cy[5], cz[5] derived from precomputed control points
   // double f = c0 + c1 * U + c2 * U2 +   c3 * U3 +    c4 * U4;
   // df =                c1   + 2*c2 * U  + 3*c3 * U2 +  4*c4 * U3;
   // ddf =                      2*c2      + 6*c3 * U  + 12*c4 * U2;
   const array<double,5> U{1,pow(u,1),pow(u,2),pow(u,3),pow(u,4)};
   // initialize 
   for(int i=0;i<2;i++){
      v[i]=0;
      dv[i]=0;
      ddv[i]=0;
   }
   // compute curve
   for (int i=0;i<5;i++){
      v[0]+=cx[i]*U[i];
      v[1]+=cy[i]*U[i];
   }
   // compute 1st derivative of curve wrt u
   for (int i=1;i<5;i++){
      dv[0]+=i*cx[i]*U[i-1];
      dv[1]+=i*cy[i]*U[i-1];
   }
   // compute 2nd derivative of curve wrt u
   for (int i=2;i<5;i++){
      ddv[0]+=i*(i-1)*cx[i]*U[i-2];
      ddv[1]+=i*(i-1)*cy[i]*U[i-2];
   }
}

// Function ptr for 1d computes, used for 1d root finding in following function
typedef void(*fPtr)(const array<double,5>& coefficients, double u, array<double,3>& c);

// Newton solver; find root; compute function defined above
double get_parametric_coordinate_at_max(const array<double,5>& cy, double u0, double tolerance, int max_iter)
{
   //const double u0=0.5;
   double u=u0;
   array<double,3> v;
   compute(cy,u,v);
   // 1st and 2nd derivative
   double f=v[1];
   double j=v[2];
   double resid=std::fabs(f);
   int iter=0;
   while (resid>tolerance && iter<=max_iter){
      // Compute increment and update newton iteration
      // u+du = inv(Jacobian) * residual
      u+=-f/j;
      compute(cy,u,v);
      f=v[1];
      j=v[2];
      resid=std::fabs(f);
      iter+=1;
   }
   return u;
}

double 
cpp_predictor_2d(const array<double,5>& cx, const array<double,5>& cy, 
             const array<double,2>& xy, double umax_width, double du)
{

   using fourth_order_bezier_polynomial::compute;
   using fourth_order_bezier_polynomial::curve;
   const double tolerance=1.0e-15;
   double x=xy[0];
   double y=xy[1];
   double xm;
   {
      // Compute x-coordinate at maximum y-value
      xm=compute(cx,umax_width);
   }
   //printf("predictor_2d umax_depth=%12.6f\n",umax_width);
   //printf("predictor_2d xm=%12.6f\n",xm);
   double s;
   double umin, umax, minlim, maxlim;
   // Identify quadrant of input spatial point (x,y)
   // Starting value for bounds on parametric coordinate
   minlim=0.0;
   maxlim=1.0;
   if (x>xm){
      umin=umax_width;
      umax=1.0;
   } else {
      umin=0;
      umax=umax_width;
   }

   {
      // Bound parametic coordinate and compute (umin, umax) of bounds
      double un=0.5*(umin+umax);
      while (umax-umin > du){
         array<double,2> rho;
         array<double,2> rho1;
         curve(cx,cy,un,rho,rho1);
         array<double,2> dr{x-rho[0],y-rho[1]};
         double dt=dot(dr,rho1);
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

      //printf("predictor_2d umax=%12.6f\n",umax);
      //printf("predictor_2d umin=%12.6f\n",umin);
      //printf("predictor_2d un=%12.6f\n",un);
   }

   {
      // Compute points on curve at c0(umin), c1(umax)
      array<double,2> c0;
      array<double,2> c1;
      curve(cx,cy,umin,c0);
      curve(cx,cy,umax,c1);
      array<double,2> delta{c0[0]-c1[0],c0[1]-c1[1]};
      array<double,2> dr{x-c0[0],y-c0[1]};
      double t=dot(delta,dr)/dot(delta,delta);
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
   //printf("predictor_2d s=%12.6f\n",s);
   return s;
}

double 
solve(const array<double,5>& cx, double u0, double x, double tolerance, int max_iter) 
{
   // Given an 'x', solve for u; 
   // u0 is starting guess
   double u=u0;
   array<double,3> curve;
   compute(cx,u,curve);
   double f=-(curve[0]-x);
   double j=curve[1];
   double resid=std::fabs(f);
   int iter=0;
   while (resid>tolerance && iter<=max_iter){
      // Compute increment and update newton iteration
      // u+du = inv(Jacobian) * residual
      double du=f/j;
      while(u+du<0 || u+du>1) du*=0.5;
      u+=du;
      compute(cx,u,curve);
      f=-(curve[0]-x);
      j=curve[1];
      resid=std::fabs(f);
      iter+=1;
   }
   return u;
}

tuple<double,double> 
closest_point(const array<double,5>& cx, const array<double,5>& cy, 
    const array<double,2>& xy, double umax_curve, double tolerance, int max_iter,
    PathDot path_dot) {
   // returns tuple<double,double>(u, distance)
   // Compute parametric coordinate defining closest point;  
   // Using parametric coordinate for closest point on curve, 
   // function also returns distance to curve.
   //
   // Initial guess must be decent otherwise cpp fails
   double du=0.025;
   double u0=cpp_predictor_2d(cx,cy,xy,umax_curve,du);

   // Initial estimate
   double u=u0;
   array<double,2> rho{0,0};
   array<double,2> rho_1{0,0};
   array<double,2> rho_11{0,0};
   // evaluate curve
   curve(cx,cy,u,rho,rho_1,rho_11);
   array<double,2> dr{xy[0]-rho[0],xy[1]-rho[1]};
   {
      // Is this point an interior point?  
      // 'rho_1' is tangent vector
      // Define normal 'n' to curve so that it points towards exterior of pool; 
      // For 2d curves, CCW direction is right handed
      // For 2d curves, CW direction is left handed;
      //   for CW curves, use a dot product which is negated 
      array<double,2> n{ rho_1[1], -rho_1[0]};
      //  dr is vector pointing from 'curve' to input point 'x,y' ;
      //  negative from dot product means that dr is pointing 
      //  towards interior
      //double dt=(n[0]*dr[0]+n[1]*dr[1]);
      double dt=path_dot(n,dr);
      //printf("INITIAL cpp: @u0=%10.6f,rho={%10.6f, %10.6f}\n",u0,rho[0],rho[1]);
      //printf("INITIAL cpp: @u0=%10.6f, dr={%10.6f, %10.6f}\n",u0,dr[0],dr[1]);
      //printf("INITIAL cpp: @u0=%10.6f,  n={%10.6f, %10.6f}\n",u0,n[0],n[1]);
      //printf("INITIAL cpp: @u0=%10.6f, dt = dot(n,dr) = %10.6f\n",u0,dt);
      // If above dot product is negative then (x,y) defines 
      // a point inside the pool.
      if (dt<=0) {
         return std::make_tuple(0.0,-1.0);
      }
   }
   // compute_ddcurve(u,rho_11);
   // Jacobian
   double j11=dot(rho_11,dr)-dot(rho_1,rho_1);
   //std::cout << "dr = " << std::scientific << std::setw(15) << dr[0] << ", " << dr[1] << std::endl;
   //std::cout << "rho_1 = " << std::scientific << std::setw(15) << rho_1[0] << ", " << rho_1[1] << std::endl;

   // Compute residual
   double r1=dot(rho_1,dr);
   double resid=std::fabs(r1);

   int iter=0;
   //std::cout << "Iteration # Solution Value  Residual" << std::endl;
   //std::cout << "----------- --------------  --------" << std::endl;
   //std::cout << std::setw(10) << iter;
   //std::cout << std::scientific << std::setw(15) << u;
   //std::cout << std::scientific << std::setw(15) << resid << std::endl;
   while (resid>tolerance && iter<=max_iter){

      // Compute increment and update newton iteration
      // u+du = inv(Jacobian) * residual
      u+=-r1/j11;

      // Update residual in preparation for next iteration 
      curve(cx,cy,u,rho,rho_1,rho_11);
      //compute_curve(u,rho,rho_1,rho_11);
      //compute_curve(u,rho);
      dr[0]=xy[0]-rho[0];
      dr[1]=xy[1]-rho[1];
      // compute_dcurve(u,rho_1);
      r1=dot(rho_1,dr);
      resid=std::fabs(r1);

      // In preparation of next Jacobian calculation
      // compute_ddcurve(u,rho_11);
      // Jacobian
      j11=dot(rho_11,dr)-dot(rho_1,rho_1);

      // Go to next iteration
      iter+=1;
      //std::cout << std::setw(10) << iter;
      //std::cout << std::scientific << std::setw(15) << u;
      //std::cout << std::scientific << std::setw(15) << resid << std::endl;
   }

   // Is this point an interior point?  
   // 'rho_1' is tangent vector
   // Define normal 'n' to curve so that it points towards exterior of pool; 
   // For 2d curves, CCW direction is right handed
   // For 2d curves, CW direction is left handed;
   //   for CW curves, use a dot product which is negated 
   array<double,2> n{ rho_1[1], -rho_1[0]};
   double distance=sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
   //  dr is vector pointing from 'curve' to input point 'x,y' ;
   //  negative from dot product means that dr is pointing 
   //  towards interior
   //double dt=(n[0]*dr[0]+n[1]*dr[1]);
   double dt=path_dot(n,dr);
   //printf("FINAL cpp: @u=%10.6f,rho={%10.6f, %10.6f}\n",u,rho[0],rho[1]);
   //printf("FINAL cpp: @u=%10.6f, dr={%10.6f, %10.6f}\n",u,dr[0],dr[1]);
   //printf("FINAL cpp: @u=%10.6f,  n={%10.6f, %10.6f}\n",u,n[0],n[1]);
   //printf("FINAL cpp: @u=%10.6f, dt = dot(n,dr) = %10.6f\n",u,dt);
   if (dt<=0) distance=-1;
   return make_tuple(u,distance);

}

tuple<double, array<double,5>, array<double,5>>
scale_coordinates_for_curve_height(const array<double,5>&xt, const array<double,5>&yt, double height)
{
   // Scale top curve control points (xt,yt) for height
   array<double,5> cy=compute_coefficients(yt);
   // Newton solve to parametric coordinate for maximum height
   const double u0=0.5;
   const int max_iter=15;
   const double tolerance=1.0e-15;
   // Computing private member variable
   double umax_height=get_parametric_coordinate_at_max(cy,u0,tolerance,max_iter);
   double ymax=std::fabs(compute(cy,umax_height));
   double ymin=0;
   // initial height
   double h0=ymax-ymin;
   double alpha=height/h0;
   // scale control points of pool for height
   // scaling private variable control points for height
   array<double,5> x,y;
   for(int i=0;i<5;i++){
      x[i]=xt[i]*alpha;
      y[i]=yt[i]*alpha;
   }
   return make_tuple(umax_height,x,y);
}

}
