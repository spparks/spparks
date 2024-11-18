/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.

   quaternion.h
   Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
   Date: March 8 2024

------------------------------------------------------------------------- */
#ifndef SPK_QUATERNION_H
#define SPK_QUATERNION_H

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <math.h>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

using std::vector;

namespace SPPARKS_NS {

namespace quaternion {

inline unsigned seed =
    std::chrono::system_clock::now().time_since_epoch().count();
inline std::default_random_engine generator(seed);

inline double norm(const vector<double> &q) {
  // if (q.size() != 4)
  //   error->all(FLERR, "vector storing quaternion must have size 4.");

  double init = 0.0;
  auto r = std::inner_product(q.cbegin(), q.cend(), q.cbegin(), init);
  return std::sqrt(r);
}

inline bool check_validity(const vector<double> &q, double epsilon = 1.0e-15) {
  double n = norm(q);
  if (std::fabs(n - 1.0) > epsilon)
    throw std::runtime_error(
        "FAILS unit quaternion test: abs(||q||-1.0)>1.0e-15");
  return true;
}

inline vector<double> vector_cross_product(const vector<double> &q,
                                           const vector<double> &r) {
  if (3 != q.size() || 3 != r.size())
    throw std::runtime_error("len(q)!=3 or len(r)!=3");
  double qx = q[0], qy = q[1], qz = q[2];
  double rx = r[0], ry = r[1], rz = r[2];
  vector<double> s{qy * rz - qz * ry, -(qx * rz - qz * rx), qx * ry - qy * rx};
  return std::move(s);
}
inline vector<double> composition(const vector<double> &a,
                                  const vector<double> &b, bool unit = true) {
  // Computes quaternion composition a b
  //   args: a quaternion array len=4.
  //   args: b quaternion array len=4.
  //   if unit=true (default) a and b must be unit quaternions
  //   returns: quaternion array len=4;

  // if (a.size() != 4 || b.size() != 4)
  //   error->all(FLERR, "vector storing quaternion must have size 4.");

  if (unit) {
    check_validity(a);
    check_validity(b);
  }
  double ar = a[0];
  double br = b[0];
  double ax = a[1], ay = a[2], az = a[3];
  double bx = b[1], by = b[2], bz = b[3];
  // scalar qr=ar br - a[1:] dot b[1:]
  double qr = ar * br - (ax * bx + ay * by + az * bz);
  vector<double> cross{ay * bz - az * by, -(ax * bz - az * bx),
                       ax * by - ay * bx};
  // 3d vector qv=ar b[1:] + br a[1:] - a[1:] cross b[1:]
  vector<double> qv{ar * bx + br * ax + cross[0], ar * by + br * ay + cross[1],
                    ar * bz + br * az + cross[2]};
  return std::move(vector<double>{qr, qv[0], qv[1], qv[2]});
}

inline vector<double> conjugate(const vector<double> &q) {
  // if (q.size() != 4)
  //   error->all(FLERR, "vector storing quaternion must have size 4.");
  check_validity(q);
  return std::move(vector<double>{q[0], -q[1], -q[2], -q[3]});
}

inline vector<double> rotate_vector(const vector<double> &q,
                                    const vector<double> &r) {
  // """
  // Rotates vector 'r' with input quaternion q.
  // args: q unit quaternion array len=4.
  // args: r array len=3.
  // return value array len=3
  // Vec(r)+2 Real(q)Vec(q) CROSS Vec(r)+2 Vec(q) CROSS (Vec(q) CROSS Vec(r))
  // """
  if (4 != q.size() || 3 != r.size())
    throw std::runtime_error("len(q)!=4 or len(r)!=3");
  check_validity(q);
  // Real and vector parts of q
  double qr = q[0];
  vector<double> qv{q[1], q[2], q[3]};
  // rot=r+2*qr*c+2*numpy.cross(qv,c)
  vector<double> c = vector_cross_product(qv, r);
  vector<double> c2 = vector_cross_product(qv, c);
  vector<double> rot{r[0] + 2.0 * qr * c[0] + 2.0 * c2[0],
                     r[1] + 2.0 * qr * c[1] + 2.0 * c2[1],
                     r[2] + 2.0 * qr * c[2] + 2.0 * c2[2]};
  return std::move(rot);
}

inline vector<double> to_rotation_matrix(const vector<double> &q) {
  //     """
  //     Converts quaternion q to 3x3 unitary matrix
  //     returns 'row major' 3x3 matrix.
  //
  //     Return value: Columns of matrix express local site orientation
  //     (per input quaternion) with respect to problem
  //     global coordinate axes.  Depending on context, it may
  //     be necessary to transpose this matrix.
  //     """
  //     check_validity(q,unit=True)

  // if (q.size() != 4)
  //   error->all(FLERR, "vector storing quaternion must have size 4.");
  const double q0 = q[0], qx = q[1], qy = q[2], qz = q[3];
  const double q00 = q0 * q0;
  const double q0x = q0 * qx;
  const double q0y = q0 * qy;
  const double q0z = q0 * qz;
  const double qxx = qx * qx;
  const double qxy = qx * qy;
  const double qxz = qx * qz;
  const double qyy = qy * qy;
  const double qyz = qy * qz;
  const double qzz = qz * qz;
  const double r11 = q00 + qxx - qyy - qzz;
  const double r22 = q00 - qxx + qyy - qzz;
  const double r33 = q00 - qxx - qyy + qzz;
  const double r12 = 2 * (qxy - q0z);
  const double r21 = 2 * (qxy + q0z);
  const double r13 = 2 * (qxz + q0y);
  const double r31 = 2 * (qxz - q0y);
  const double r23 = 2 * (qyz - q0x);
  const double r32 = 2 * (qyz + q0x);
  vector<double> m{r11, r12, r13, r21, r22, r23, r31, r32, r33};
  return std::move(m);
}

inline double trace_of_rotation_matrix(const vector<double> &q) {
  //     """
  //     Converts quaternion q to 3x3 unitary matrix m
  //     returns trace(m)
  //     """
  //     check_validity(q,unit=True)

  // if (q.size() != 4)
  //   error->all(FLERR, "vector storing quaternion must have size 4.");
  check_validity(q);
  const double q0 = q[0], qx = q[1], qy = q[2], qz = q[3];
  const double q00 = q0 * q0;
  const double qxx = qx * qx;
  const double qyy = qy * qy;
  const double qzz = qz * qz;
  const double r11 = q00 + qxx - qyy - qzz;
  const double r22 = q00 - qxx + qyy - qzz;
  const double r33 = q00 - qxx - qyy + qzz;
  return r11 + r22 + r33;
}

inline vector<double> to_rodrigues_vector(const vector<double> &q,
                                          double epsilon = 1.0e-4) {
  // Converts quaternion q to 3-component Rodrigues vector.
  // see also https://doi.org/10.1107/S0108767391006864
  // If q=(cos(theta/2),sin(theta/2)*[nx,ny,nz])
  //    where n = nx,ny,nz is a unit vector, then
  // the Rodrigues vector  of q is:
  // d := (1/q_0)(q_1; q_2; q_3) = tan(theta/2)*n
  //    where tan(theta/2) = sin(theta/2)/cos(theta/2)
  //    and n is the same unit vector defined above.
  // Returns a vector of length 3.

  // confirm vector length
  if (4 != q.size())
    throw std::runtime_error("len(q)!=4");

  // assign quaternion
  const double q0 = q[0], qx = q[1], qy = q[2], qz = q[3];

  // find norm via angle (see Heinz 1991 eqns 3 & 7)
  const double rho = acos(q0); // theta*2
  const double rtan = tan(rho);
  const double norm = (q0 * rtan);

  // calculate and return Rodrigues vector
  if (norm < epsilon) {
    // if angle is nearing 0
    return std::move(vector<double>{0, 0, 0});
  } else {
    // convert quaternion components to Rodrigues components
    const double rx = qx / norm;
    const double ry = qy / norm;
    const double rz = qz / norm;

    return std::move(vector<double>{rx, ry, rz});
  };
}

inline vector<double>
generate_random_unit_quaternions(std::size_t n, double epsilon = 1.0e-15) {
  // Computes n random unit quaternions
  // returns vector length = 4 * n
  // ith quaternion uq[4*i:4*i+4] for i=0,1,2,...n
  vector<double> uq(4 * n);

  // Gaussian distribution
  double mean(0.0), std(1.0);
  std::normal_distribution<double> rng(mean, std);
  for (std::size_t i = 0; i < n; i++) {
    double mag = 0.0;
    vector<double> qi{0, 0, 0, 0};
    // 4 random numbers from normal (Gaussian) distribution
    // Make sure to generate non-trivial (non-zero) quaternion
    while (mag < epsilon) {
      double sq = 0.0;
      for (std::size_t j = 0; j < 4; j++) {
        double dj = rng(generator);
        qi[j] = dj;
        sq += dj * dj;
      }
      mag = sqrt(sq);
    }

    // Normalize random quaternion
    double d = 1.0 / mag;
    for (std::size_t j = 0; j < 4; j++) {
      qi[j] *= d;
    }

    // Equivalent to above normalization loop.
    // Normalize using less obvious c++ lambda function
    //   Requires #include<algorithm>
    // std::transform(qi.begin(), qi.end(), qi.begin(),
    //                [&d](auto &c) { return d * c; });

    // Copy normalized quaternion to array of unit quaternions
    check_validity(qi, epsilon);
    for (std::size_t j = 0; j < 4; j++)
      uq[4 * i + j] = qi[j];
  }
  return std::move(uq);
}

inline double
get_cosine_of_minumum_angle_between_q_and_u(const vector<double> &q,
                                            const vector<double> &u) {
  /**
  arg q: array shape=(4,) unit quaterion
  arg u: array shape=(3,) unit vector

  Measure/compute cosine of angle between input unit quaterion 'q'
  representing a cube and the direction 'u'.  Input vector makes 3 different
  angles with local cubic coordinate axes attached to cube.  The cubic
  coordinate axes are normal to each of the faces.  Calculation of angle
  makes no distinction between positive or negative local coordinate axes of
  cube.   Intended application is to determine angle between u and line taken
  normal to face of cube where u is temperature gradient.

  returns cos(angle)
  */
  if (4 != q.size() || 3 != u.size())
    throw std::runtime_error("len(q)!=4 or len(u)!=3");

  // Take a copy of 'u' and normalize
  vector<double> n{u[0], u[1], u[2]};
  {
    double d = norm(n);
    double b = 1 / d;
    n[0] *= b;
    n[1] *= b;
    n[2] *= b;
  }
  vector<double> m = to_rotation_matrix(q);
  // Dot columns of 'm' with 'n'; which is a matrix
  // multiplication of transpose(m) X n
  double max = 0.0;
  // Because 'm' is an array in 'row-major' form -- need
  // for matrix multiply 'transpose(m)' times 'n'.
  for (size_t col = 0; col < 3; col++) {
    double dot = 0.0;
    for (size_t row = 0; row < 3; row++) {
      size_t ptr = col + 3 * row;
      dot += m[ptr] * n[row];
    }
    // Take absolute value
    dot = std::fabs(dot);
    if (dot > max)
      max = dot;
  }
  return max;
}

} // namespace quaternion

} // namespace SPPARKS_NS
#endif
