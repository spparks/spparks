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

#ifndef SPK_MATH_EXTRA_H
#define SPK_MATH_EXTRA_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "error.h"

namespace MathExtra {

  // 3 vector operations

  inline void norm3(double *v);
  inline void normalize3(const double *v, double *ans);
  inline void snormalize3(const double, const double *v, double *ans);
  inline void negate3(double *v);
  inline void scale3(double s, double *v);
  inline void add3(const double *v1, const double *v2, double *ans);
  inline void sub3(const double *v1, const double *v2, double *ans);
  inline double len3(const double *v);
  inline double lensq3(const double *v);
  inline double dot3(const double *v1, const double *v2);
  inline void cross3(const double *v1, const double *v2, double *ans);

  // 3x3 matrix operations

  inline double det3(const double mat[3][3]);
  inline void diag_times3(const double *diagonal, const double mat[3][3],
                          double ans[3][3]);
  inline void plus3(const double m[3][3], const double m2[3][3],
                    double ans[3][3]);
  inline void times3(const double m[3][3], const double m2[3][3],
                     double ans[3][3]);
  inline void transpose_times3(const double mat1[3][3], 
                               const double mat2[3][3],
                               double ans[3][3]);
  inline void times3_transpose(const double mat1[3][3], 
			       const double mat2[3][3],
			       double ans[3][3]);
  inline void invert3(const double mat[3][3], double ans[3][3]);
  inline void matvec(const double mat[3][3], const double*vec, double *ans);
  inline void matvec(const double *ex, const double *ey, const double *ez,
		     const double *vec, double *ans);
  inline void transpose_matvec(const double mat[3][3], const double*vec,
			       double *ans);
  inline void transpose_matvec(const double *ex, const double *ey, 
			       const double *ez, const double *v,
			       double *ans);
  inline void transpose_diag3(const double mat[3][3], const double*vec,
			      double ans[3][3]);
  inline void vecmat(const double *v, const double m[3][3], double *ans);
  inline void scalar_times3(const double f, double m[3][3]); 

  void write3(const double mat[3][3]);
  int mldivide3(const double mat[3][3], const double *vec, double *ans);
  int jacobi(double matrix[3][3], double *evalues, double evectors[3][3]);
  void rotate(double matrix[3][3], int i, int j, int k, int l,
	      double s, double tau);
  void richardson(double *q, double *m, double *w, double *moments, double dtq);

  // shape matrix operations
  // upper-triangular 3x3 matrix stored in Voigt notation as 6-vector

  inline void multiply_shape_shape(const double *one, const double *two,
                                   double *ans);

  // quaternion operations

  inline void qnormalize(double *q);
  inline void qconjugate(double *q, double *qc);
  inline void vecquat(double *a, double *b, double *c);
  inline void quatvec(double *a, double *b, double *c);
  inline void quatquat(double *a, double *b, double *c);
  inline void invquatvec(double *a, double *b, double *c);
  inline void axisangle_to_quat(const double *v, const double angle,
                                double *quat);

  void angmom_to_omega(double *m, double *ex, double *ey, double *ez,
		       double *idiag, double *w);
  void omega_to_angmom(double *w, double *ex, double *ey, double *ez,
		       double *idiag, double *m);
  void mq_to_omega(double *m, double *q, double *moments, double *w);
  void exyz_to_q(double *ex, double *ey, double *ez, double *q);
  void q_to_exyz(double *q, double *ex, double *ey, double *ez);
  void quat_to_mat(const double *quat, double mat[3][3]);
  void quat_to_mat_trans(const double *quat, double mat[3][3]);

  // rotation operations
  
  inline void rotation_generator_x(const double m[3][3], double ans[3][3]);
  inline void rotation_generator_y(const double m[3][3], double ans[3][3]);
  inline void rotation_generator_z(const double m[3][3], double ans[3][3]);

  // moment of inertia operations

  void inertia_ellipsoid(double *shape, double *quat, double mass,
			 double *inertia);
  void inertia_line(double length, double theta, double mass,
		    double *inertia);
  void inertia_triangle(double *v0, double *v1, double *v2, 
			double mass, double *inertia);
  void inertia_triangle(double *idiag, double *quat, double mass, 
			double *inertia); 
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

void MathExtra::norm3(double *v)
{
  double scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

void MathExtra::normalize3(const double *v, double *ans)
{
  double scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

void MathExtra::snormalize3(const double length, const double *v, double *ans)
{
  double scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   negate vector v
------------------------------------------------------------------------- */

void MathExtra::negate3(double *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s
------------------------------------------------------------------------- */

void MathExtra::scale3(double s, double *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

void MathExtra::add3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

void MathExtra::sub3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

double MathExtra::len3(const double *v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

double MathExtra::lensq3(const double *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

double MathExtra::dot3(const double *v1, const double *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

void MathExtra::cross3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[1]*v2[2]-v1[2]*v2[1];
  ans[1] = v1[2]*v2[0]-v1[0]*v2[2];
  ans[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

double MathExtra::det3(const double m[3][3])
{
  double ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] - 
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] + 
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

void MathExtra::diag_times3(const double *d, const double m[3][3],
			    double ans[3][3])
{
  ans[0][0] = d[0]*m[0][0];
  ans[0][1] = d[0]*m[0][1];
  ans[0][2] = d[0]*m[0][2];
  ans[1][0] = d[1]*m[1][0];
  ans[1][1] = d[1]*m[1][1];
  ans[1][2] = d[1]*m[1][2];
  ans[2][0] = d[2]*m[2][0];
  ans[2][1] = d[2]*m[2][1];
  ans[2][2] = d[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

void MathExtra::plus3(const double m[3][3], const double m2[3][3],
		      double ans[3][3])
{
  ans[0][0] = m[0][0]+m2[0][0];
  ans[0][1] = m[0][1]+m2[0][1];
  ans[0][2] = m[0][2]+m2[0][2];
  ans[1][0] = m[1][0]+m2[1][0];
  ans[1][1] = m[1][1]+m2[1][1];
  ans[1][2] = m[1][2]+m2[1][2];
  ans[2][0] = m[2][0]+m2[2][0];
  ans[2][1] = m[2][1]+m2[2][1];
  ans[2][2] = m[2][2]+m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

void MathExtra::times3(const double m[3][3], const double m2[3][3],
                       double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[1][0] + m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[0][1]*m2[1][1] + m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[0][1]*m2[1][2] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[1][0] + m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1] + m[1][1]*m2[1][1] + m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2] + m[1][1]*m2[1][2] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1] + m[2][1]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2] + m[2][1]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

void MathExtra::transpose_times3(const double m[3][3], const double m2[3][3],
                                 double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[1][0]*m2[1][0] + m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[1][0]*m2[1][1] + m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[1][0]*m2[1][2] + m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0] + m[1][1]*m2[1][0] + m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1] + m[1][1]*m2[1][1] + m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2] + m[1][1]*m2[1][2] + m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0] + m[1][2]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1] + m[1][2]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2] + m[1][2]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

void MathExtra::times3_transpose(const double m[3][3], const double m2[3][3],
                                 double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[0][1] + m[0][2]*m2[0][2];
  ans[0][1] = m[0][0]*m2[1][0] + m[0][1]*m2[1][1] + m[0][2]*m2[1][2];
  ans[0][2] = m[0][0]*m2[2][0] + m[0][1]*m2[2][1] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[0][1] + m[1][2]*m2[0][2];
  ans[1][1] = m[1][0]*m2[1][0] + m[1][1]*m2[1][1] + m[1][2]*m2[1][2];
  ans[1][2] = m[1][0]*m2[2][0] + m[1][1]*m2[2][1] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[0][1] + m[2][2]*m2[0][2];
  ans[2][1] = m[2][0]*m2[1][0] + m[2][1]*m2[1][1] + m[2][2]*m2[1][2];
  ans[2][2] = m[2][0]*m2[2][0] + m[2][1]*m2[2][1] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT checks for singular or badly scaled matrix
------------------------------------------------------------------------- */

void MathExtra::invert3(const double m[3][3], double ans[3][3])
{
  double den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
  den += -m[1][0]*m[0][1]*m[2][2]+m[1][0]*m[0][2]*m[2][1];
  den += m[2][0]*m[0][1]*m[1][2]-m[2][0]*m[0][2]*m[1][1];

  ans[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) / den;
  ans[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1]) / den;
  ans[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) / den;
  ans[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0]) / den;
  ans[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) / den;
  ans[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0]) / den;
  ans[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) / den;
  ans[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0]) / den;
  ans[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

void MathExtra::matvec(const double m[3][3], const double *v, double *ans) 
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

void MathExtra::matvec(const double *ex, const double *ey, const double *ez,
		       const double *v, double *ans) 
{
  ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
  ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
  ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

void MathExtra::transpose_matvec(const double m[3][3], const double *v,
				 double *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

void MathExtra::transpose_matvec(const double *ex, const double *ey, 
				 const double *ez, const double *v,
				 double *ans)
{
  ans[0] = ex[0]*v[0] + ex[1]*v[1] + ex[2]*v[2];
  ans[1] = ey[0]*v[0] + ey[1]*v[1] + ey[2]*v[2];
  ans[2] = ez[0]*v[0] + ez[1]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

void MathExtra::transpose_diag3(const double m[3][3], const double *d, 
				double ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[1][0]*d[1];
  ans[0][2] = m[2][0]*d[2];
  ans[1][0] = m[0][1]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[2][1]*d[2];
  ans[2][0] = m[0][2]*d[0];
  ans[2][1] = m[1][2]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

void MathExtra::vecmat(const double *v, const double m[3][3], double *ans)
{
  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

inline void MathExtra::scalar_times3(const double f, double m[3][3]) 
{
  m[0][0] *= f; m[0][1] *= f; m[0][2] *= f;
  m[1][0] *= f; m[1][1] *= f; m[1][2] *= f;
  m[2][0] *= f; m[2][1] *= f; m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   multiply 2 shape matrices
   upper-triangular 3x3, stored as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void MathExtra::multiply_shape_shape(const double *one, const double *two,
                                     double *ans)
{
  ans[0] = one[0]*two[0];
  ans[1] = one[1]*two[1];
  ans[2] = one[2]*two[2];
  ans[3] = one[1]*two[3] + one[3]*two[2];
  ans[4] = one[0]*two[4] + one[5]*two[3] + one[4]*two[2];
  ans[5] = one[0]*two[5] + one[5]*two[1];
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

void MathExtra::qnormalize(double *q)
{
  double norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

void MathExtra::qconjugate(double *q, double *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   vector-quaternion multiply: c = a*b, where a = (0,a)
------------------------------------------------------------------------- */

void MathExtra::vecquat(double *a, double *b, double *c)
{
  c[0] = -a[0]*b[1] - a[1]*b[2] - a[2]*b[3];
  c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
  c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
  c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
}

/* ----------------------------------------------------------------------
   quaternion-vector multiply: c = a*b, where b = (0,b)
------------------------------------------------------------------------- */

void MathExtra::quatvec(double *a, double *b, double *c)
{
  c[0] = -a[1]*b[0] - a[2]*b[1] - a[3]*b[2];
  c[1] = a[0]*b[0] + a[2]*b[2] - a[3]*b[1];
  c[2] = a[0]*b[1] + a[3]*b[0] - a[1]*b[2];
  c[3] = a[0]*b[2] + a[1]*b[1] - a[2]*b[0];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */

void MathExtra::quatquat(double *a, double *b, double *c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  c[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
}

/* ----------------------------------------------------------------------
   quaternion multiply: c = inv(a)*b
   a is a quaternion
   b is a four component vector
   c is a three component vector
------------------------------------------------------------------------- */

void MathExtra::invquatvec(double *a, double *b, double *c)
{ 
  c[0] = -a[1]*b[0] + a[0]*b[1] + a[3]*b[2] - a[2]*b[3];
  c[1] = -a[2]*b[0] - a[3]*b[1] + a[0]*b[2] + a[1]*b[3];
  c[2] = -a[3]*b[0] + a[2]*b[1] - a[1]*b[2] + a[0]*b[3];
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

void MathExtra::axisangle_to_quat(const double *v, const double angle,
				  double *quat)
{
  double halfa = 0.5*angle;
  double sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about x to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_x(const double m[3][3], double ans[3][3])
{
  ans[0][0] = 0;
  ans[0][1] = -m[0][2];
  ans[0][2] = m[0][1];
  ans[1][0] = 0;
  ans[1][1] = -m[1][2];
  ans[1][2] = m[1][1];
  ans[2][0] = 0;
  ans[2][1] = -m[2][2];
  ans[2][2] = m[2][1];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about y to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_y(const double m[3][3], double ans[3][3])
{
  ans[0][0] = m[0][2];
  ans[0][1] = 0;
  ans[0][2] = -m[0][0];
  ans[1][0] = m[1][2];
  ans[1][1] = 0;
  ans[1][2] = -m[1][0];
  ans[2][0] = m[2][2];
  ans[2][1] = 0;
  ans[2][2] = -m[2][0];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about z to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_z(const double m[3][3], double ans[3][3])
{
  ans[0][0] = -m[0][1];
  ans[0][1] = m[0][0];
  ans[0][2] = 0;
  ans[1][0] = -m[1][1];
  ans[1][1] = m[1][0];
  ans[1][2] = 0;
  ans[2][0] = -m[2][1];
  ans[2][1] = m[2][0];
  ans[2][2] = 0;
}

#endif
