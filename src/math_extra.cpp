/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "math_extra.h"

#define MAXJACOBI 50

namespace MathExtra {

/* ----------------------------------------------------------------------
   output a matrix
------------------------------------------------------------------------- */

void write3(const double mat[3][3])
{
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) printf("%g ",mat[i][j]);
    printf("\n");
  }
}

/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
------------------------------------------------------------------------- */

int mldivide3(const double m[3][3], const double *v, double *ans)
{
  // create augmented matrix for pivoting

  double aug[3][4];
  for (unsigned i = 0; i < 3; i++) {
    aug[i][3] = v[i];
    for (unsigned j = 0; j < 3; j++) aug[i][j] = m[i][j];
  }

  for (unsigned i = 0; i < 2; i++) {
    unsigned p = i;
    for (unsigned j = i+1; j < 3; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memcpy(aug[i],aug[j],4*sizeof(double));
        memcpy(aug[j],tempv,4*sizeof(double));
      }
    }

    while (aug[p][i] == 0.0 && p < 3) p++;

    if (p == 3) return 1;
    else
      if (p != i) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memcpy(aug[i],aug[p],4*sizeof(double));
        memcpy(aug[p],tempv,4*sizeof(double));
      }

    for (unsigned j = i+1; j < 3; j++) {
      double m = aug[j][i]/aug[i][i];
      for (unsigned k=i+1; k<4; k++) aug[j][k]-=m*aug[i][k];
    }
  }

  if (aug[2][2] == 0.0) return 1;
  
  // back substitution

  ans[2] = aug[2][3]/aug[2][2];
  for (int i = 1; i >= 0; i--) {
    double sumax = 0.0;
    for (unsigned j = i+1; j < 3; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][3]-sumax) / aug[i][i];
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute evalues and evectors of 3x3 real symmetric matrix
   based on Jacobi rotations
   adapted from Numerical Recipes jacobi() function
------------------------------------------------------------------------- */

int jacobi(double matrix[3][3], double *evalues, double evectors[3][3])
{
  int i,j,k;
  double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) evectors[i][j] = 0.0;
    evectors[i][i] = 1.0;
  }
  for (i = 0; i < 3; i++) {
    b[i] = evalues[i] = matrix[i][i];
    z[i] = 0.0;
  }
  
  for (int iter = 1; iter <= MAXJACOBI; iter++) {
    sm = 0.0;
    for (i = 0; i < 2; i++)
      for (j = i+1; j < 3; j++)
	sm += fabs(matrix[i][j]);
    if (sm == 0.0) return 0;
    
    if (iter < 4) tresh = 0.2*sm/(3*3);
    else tresh = 0.0;
    
    for (i = 0; i < 2; i++) {
      for (j = i+1; j < 3; j++) {
	g = 100.0*fabs(matrix[i][j]);
	if (iter > 4 && fabs(evalues[i])+g == fabs(evalues[i])
	    && fabs(evalues[j])+g == fabs(evalues[j]))
	  matrix[i][j] = 0.0;
	else if (fabs(matrix[i][j]) > tresh) {
	  h = evalues[j]-evalues[i];
	  if (fabs(h)+g == fabs(h)) t = (matrix[i][j])/h;
	  else {
	    theta = 0.5*h/(matrix[i][j]);
	    t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt(1.0+t*t);
	  s = t*c;
	  tau = s/(1.0+c);
	  h = t*matrix[i][j];
	  z[i] -= h;
	  z[j] += h;
	  evalues[i] -= h;
	  evalues[j] += h;
	  matrix[i][j] = 0.0;
	  for (k = 0; k < i; k++) rotate(matrix,k,i,k,j,s,tau);
	  for (k = i+1; k < j; k++) rotate(matrix,i,k,k,j,s,tau);
	  for (k = j+1; k < 3; k++) rotate(matrix,i,k,j,k,s,tau);
	  for (k = 0; k < 3; k++) rotate(evectors,k,i,k,j,s,tau);
	}
      }
    }
    
    for (i = 0; i < 3; i++) {
      evalues[i] = b[i] += z[i];
      z[i] = 0.0;
    }
  }
  return 1;
}

/* ----------------------------------------------------------------------
   perform a single Jacobi rotation
------------------------------------------------------------------------- */

void rotate(double matrix[3][3], int i, int j, int k, int l,
	    double s, double tau)
{
  double g = matrix[i][j];
  double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}

/* ----------------------------------------------------------------------
   Richardson iteration to update quaternion from angular momentum
   return new normalized quaternion q
   also returns updated omega at 1/2 step
------------------------------------------------------------------------- */

void richardson(double *q, double *m, double *w, double *moments, double dtq)
{
  // full update from dq/dt = 1/2 w q

  double wq[4];
  MathExtra::vecquat(w,q,wq);

  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtra::qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);

  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  MathExtra::mq_to_omega(m,qhalf,moments,w);
  MathExtra::vecquat(w,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtra::qnormalize(qhalf);

  // corrected Richardson update

  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];
  MathExtra::qnormalize(q);
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum, both in space frame
   only know Idiag so need to do M = Iw in body frame
   ex,ey,ez are column vectors of rotation matrix P
   Mbody = P_transpose Mspace
   wbody = Mbody / Idiag
   wspace = P wbody
   set wbody component to 0.0 if inertia component is 0.0
     otherwise body can spin easily around that axis
------------------------------------------------------------------------- */

void angmom_to_omega(double *m, double *ex, double *ey, double *ez,
		     double *idiag, double *w)
{
  double wbody[3];

  if (idiag[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / idiag[0];
  if (idiag[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / idiag[1];
  if (idiag[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / idiag[2];

  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in body frame
   project space-frame angular momentum onto body axes
     and divide by principal moments
------------------------------------------------------------------------- */

void mq_to_omega(double *m, double *q, double *moments, double *w)
{
  double wbody[3];
  double rot[3][3];

  MathExtra::quat_to_mat(q,rot);
  MathExtra::transpose_matvec(rot,m,wbody);
  if (moments[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] /= moments[0];
  if (moments[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] /= moments[1];
  if (moments[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] /= moments[2];
  MathExtra::matvec(rot,wbody,w);
}

/* ----------------------------------------------------------------------
   compute angular momentum from omega, both in space frame
   only know Idiag so need to do M = Iw in body frame
   ex,ey,ez are column vectors of rotation matrix P
   wbody = P_transpose wspace
   Mbody = Idiag wbody
   Mspace = P Mbody
------------------------------------------------------------------------- */

void omega_to_angmom(double *w, double *ex, double *ey, double *ez,
		     double *idiag, double *m)
{
  double mbody[3];

  mbody[0] = (w[0]*ex[0] + w[1]*ex[1] + w[2]*ex[2]) * idiag[0];
  mbody[1] = (w[0]*ey[0] + w[1]*ey[1] + w[2]*ey[2]) * idiag[1];
  mbody[2] = (w[0]*ez[0] + w[1]*ez[1] + w[2]*ez[2]) * idiag[2];

  m[0] = mbody[0]*ex[0] + mbody[1]*ey[0] + mbody[2]*ez[0];
  m[1] = mbody[0]*ex[1] + mbody[1]*ey[1] + mbody[2]*ez[1];
  m[2] = mbody[0]*ex[2] + mbody[1]*ey[2] + mbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   create unit quaternion from space-frame ex,ey,ez
   ex,ey,ez are columns of a rotation matrix
------------------------------------------------------------------------- */

void exyz_to_q(double *ex, double *ey, double *ez, double *q)
{
  // squares of quaternion components
  
  double q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0);
  double q1sq = q0sq - 0.5 * (ey[1] + ez[2]);
  double q2sq = q0sq - 0.5 * (ex[0] + ez[2]);
  double q3sq = q0sq - 0.5 * (ex[0] + ey[1]);
  
  // some component must be greater than 1/4 since they sum to 1
  // compute other components from it
  
  if (q0sq >= 0.25) {
    q[0] = sqrt(q0sq);
    q[1] = (ey[2] - ez[1]) / (4.0*q[0]);
    q[2] = (ez[0] - ex[2]) / (4.0*q[0]);
    q[3] = (ex[1] - ey[0]) / (4.0*q[0]);
  } else if (q1sq >= 0.25) {
    q[1] = sqrt(q1sq);
    q[0] = (ey[2] - ez[1]) / (4.0*q[1]);
    q[2] = (ey[0] + ex[1]) / (4.0*q[1]);
    q[3] = (ex[2] + ez[0]) / (4.0*q[1]);
  } else if (q2sq >= 0.25) {
    q[2] = sqrt(q2sq);
    q[0] = (ez[0] - ex[2]) / (4.0*q[2]);
    q[1] = (ey[0] + ex[1]) / (4.0*q[2]);
    q[3] = (ez[1] + ey[2]) / (4.0*q[2]);
  } else if (q3sq >= 0.25) {
    q[3] = sqrt(q3sq);
    q[0] = (ex[1] - ey[0]) / (4.0*q[3]);
    q[1] = (ez[0] + ex[2]) / (4.0*q[3]);
    q[2] = (ez[1] + ey[2]) / (4.0*q[3]);
  }

  qnormalize(q);
}

/* ----------------------------------------------------------------------
   compute space-frame ex,ey,ez from current quaternion q
   ex,ey,ez = space-frame coords of 1st,2nd,3rd principal axis
   operation is ex = q' d q = Q d, where d is (1,0,0) = 1st axis in body frame
------------------------------------------------------------------------- */

void q_to_exyz(double *q, double *ex, double *ey, double *ez)
{
  ex[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  ex[1] = 2.0 * (q[1]*q[2] + q[0]*q[3]);
  ex[2] = 2.0 * (q[1]*q[3] - q[0]*q[2]);
  
  ey[0] = 2.0 * (q[1]*q[2] - q[0]*q[3]);
  ey[1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  ey[2] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
  
  ez[0] = 2.0 * (q[1]*q[3] + q[0]*q[2]);
  ez[1] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
  ez[2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */

void quat_to_mat(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;
	
  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion conjugate
   quat = [w i j k]
------------------------------------------------------------------------- */

void quat_to_mat_trans(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[1][0] = twoij-twokw;
  mat[2][0] = twojw+twoik;

  mat[0][1] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[2][1] = twojk-twoiw;
	
  mat[0][2] = twoik-twojw;
  mat[1][2] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

/* ----------------------------------------------------------------------
   compute space-frame inertia tensor of an ellipsoid
   radii = 3 radii of ellipsoid
   quat = orientiation quaternion of ellipsoid
   return symmetric inertia tensor as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void inertia_ellipsoid(double *radii, double *quat, double mass,
		       double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];
  double idiag[3];

  quat_to_mat(quat,p);
  quat_to_mat_trans(quat,ptrans);
  idiag[0] = 0.2*mass * (radii[1]*radii[1] + radii[2]*radii[2]);
  idiag[1] = 0.2*mass * (radii[0]*radii[0] + radii[2]*radii[2]);
  idiag[2] = 0.2*mass * (radii[0]*radii[0] + radii[1]*radii[1]);
  diag_times3(idiag,ptrans,itemp);
  times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}

/* ----------------------------------------------------------------------
   compute space-frame inertia tensor of a line segment in 2d
   length = length of line
   theta = orientiation of line
   return symmetric inertia tensor as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void inertia_line(double length, double theta, double mass, double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];
  double q[4],idiag[3];

  q[0] = cos(0.5*theta);
  q[1] = q[2] = 0.0;
  q[3] = sin(0.5*theta);
  MathExtra::quat_to_mat(q,p);
  MathExtra::quat_to_mat_trans(q,ptrans);
  idiag[0] = 0.0;
  idiag[1] = 1.0/12.0 * mass * length*length;
  idiag[2] = 1.0/12.0 * mass * length*length;
  MathExtra::diag_times3(idiag,ptrans,itemp);
  MathExtra::times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}

/* ----------------------------------------------------------------------
   compute space-frame inertia tensor of a triangle
   v0,v1,v2 = 3 vertices of triangle
   from http://en.wikipedia.org/wiki/Inertia_tensor_of_triangle
   inertia tensor = a/24 (v0^2 + v1^2 + v2^2 + (v0+v1+v2)^2) I - a Vt S V
   a = 2*area of tri = |(v1-v0) x (v2-v0)|
   I = 3x3 identity matrix
   V = 3x3 matrix with v0,v1,v2 as rows
   Vt = 3x3 matrix with v0,v1,v2 as columns
   S = 1/24 [2 1 1]
            [1 2 1]
            [1 1 2]
   return symmetric inertia tensor as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void inertia_triangle(double *v0, double *v1, double *v2,
		      double mass, double *inertia)
{
  double s[3][3] = {{2.0, 1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 2.0}};
  double v[3][3],sv[3][3],vtsv[3][3];
  double vvv[3],v1mv0[3],v2mv0[3],normal[3];

  v[0][0] = v0[0]; v[0][1] = v0[1]; v[0][2] = v0[2];
  v[1][0] = v1[0]; v[1][1] = v1[1]; v[1][2] = v1[2];
  v[2][0] = v2[0]; v[2][1] = v2[1]; v[2][2] = v2[2];

  times3(s,v,sv);
  transpose_times3(v,sv,vtsv);

  double sum = lensq3(v0) + lensq3(v1) + lensq3(v2);
  vvv[0] = v0[0] + v1[0] + v2[0];
  vvv[1] = v0[1] + v1[1] + v2[1];
  vvv[2] = v0[2] + v1[2] + v2[2];
  sum += lensq3(vvv);

  sub3(v1,v0,v1mv0);
  sub3(v2,v0,v2mv0);
  cross3(v1mv0,v2mv0,normal);
  double a = len3(normal);
  double inv24 = mass/24.0;

  inertia[0] = inv24*a*(sum-vtsv[0][0]);
  inertia[1] = inv24*a*(sum-vtsv[1][1]);
  inertia[2] = inv24*a*(sum-vtsv[2][2]);
  inertia[3] = -inv24*a*vtsv[1][2];
  inertia[4] = -inv24*a*vtsv[0][2];
  inertia[5] = -inv24*a*vtsv[0][1];
}

/* ----------------------------------------------------------------------
   compute space-frame inertia tensor of a triangle
   idiag = previously computed diagonal inertia tensor
   quat = orientiation quaternion of triangle
   return symmetric inertia tensor as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void inertia_triangle(double *idiag, double *quat, double mass,
		      double *inertia)
{
  double p[3][3],ptrans[3][3],itemp[3][3],tensor[3][3];

  quat_to_mat(quat,p);
  quat_to_mat_trans(quat,ptrans);
  diag_times3(idiag,ptrans,itemp);
  times3(p,itemp,tensor);
  inertia[0] = tensor[0][0];
  inertia[1] = tensor[1][1];
  inertia[2] = tensor[2][2];
  inertia[3] = tensor[1][2];
  inertia[4] = tensor[0][2];
  inertia[5] = tensor[0][1];
}

/* ---------------------------------------------------------------------- */

}
