/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.

   cubic_symmetries.h
   Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
   Date: March 8 2024

------------------------------------------------------------------------- */
#ifndef SPK_CUBIC_SYMMETRIES_H
#define SPK_CUBIC_SYMMETRIES_H

#include "quaternion.h"
#include <algorithm>
#include <cstdio>
#include <math.h>
#include <vector>
using std::vector;

namespace SPPARKS_NS {

namespace CUBIC {

const double sqrt2 = std::sqrt(2.0);

inline vector<double> get_symmetries() {
  // Cubic rotational symmetries defined using 'axis-angle' unit
  // quaternion representation.
  // Rotation angle theta, axis n.
  // q=(cos(theta/2),sin(theta/2)*[nx,ny,nz])
  //    where n = nx,ny,nz is a unit vector and is referred to as 'spindle'.
  //
  // # Constants used to define unit quaternions
  const double b = 1.0 / 2;
  const double c = 1 / std::sqrt(2.0);

  // # Identity (no rotation)
  const vector<double> q1{1, 0, 0, 0};

  // #####################################################
  // ### All spindles run through center of cube
  // ##
  // # x-axis spindles
  // # +-90 degrees, 180 degress
  const vector<double> q2{c, c, 0, 0};
  const vector<double> q3{c, -c, 0, 0};
  const vector<double> q4{0, 1, 0, 0};
  // ##
  // #y - axis spindles
  // #+ - 90 degrees, 180 degress
  const vector<double> q5{c, 0, c, 0};
  const vector<double> q6{c, 0, -c, 0};
  const vector<double> q7{0, 0, 1, 0};
  // ##
  // #z - axis spindles
  // #+ - 90 degrees, 180 degress
  const vector<double> q8{c, 0, 0, -c};
  const vector<double> q9{c, 0, 0, c};
  const vector<double> q10{0, 0, 0, 1};

  // ##
  // #Diagonal spindles between opposite edges.
  // #Cube has 12 edges with 6 pairs of opposite edges.
  // # 180 degrees about spindle.
  // #xz - plane
  const vector<double> q11{0, c, 0, c};
  const vector<double> q12{0, -c, 0, c};
  // #xy - plane
  const vector<double> q13{0, c, c, 0};
  const vector<double> q14{0, -c, c, 0};
  // #yz - plane
  const vector<double> q15{0, 0, c, c};
  const vector<double> q16{0, 0, -c, c};

  // ##
  // #Main diagonal spindles between opposite corners.
  // #Cube has 8 corners with 4 pairs of opposite corners.
  // #+ - 120 degrees about spindle.
  // #b = cos(2 * pi / 3 / 2) = cos(pi / 3) = 1 / 2
  // #n = sin(2 * pi / 3 / 2) * (1 / sqrt(3)) = 1 / 2
  // #+ 120
  const vector<double> q17{b, b, b, b};
  const vector<double> q18{b, -b, b, b};
  const vector<double> q19{b, -b, b, -b};
  const vector<double> q20{b, b, b, -b};
  // #- 120
  const vector<double> q21{b, b, -b, b};
  const vector<double> q22{b, -b, -b, b};
  const vector<double> q23{b, -b, -b, -b};
  const vector<double> q24{b, b, -b, -b};

  vector<vector<double>> symmetries(24);
  symmetries[0] = q1;
  symmetries[1] = q2;
  symmetries[2] = q3;
  symmetries[3] = q4;
  symmetries[4] = q5;
  symmetries[5] = q6;
  symmetries[6] = q7;
  symmetries[7] = q8;
  symmetries[8] = q9;
  symmetries[9] = q10;
  symmetries[10] = q11;
  symmetries[11] = q12;
  symmetries[12] = q13;
  symmetries[13] = q14;
  symmetries[14] = q15;
  symmetries[15] = q16;
  symmetries[16] = q17;
  symmetries[17] = q18;
  symmetries[18] = q19;
  symmetries[19] = q20;
  symmetries[20] = q21;
  symmetries[21] = q22;
  symmetries[22] = q23;
  symmetries[23] = q24;
  // Copy symmetries to contiguous vector
  vector<double> uq(24 * 4);
  for (std::size_t i = 0; i < 24; i++) {
    vector<double> q = symmetries[i];
    for (std::size_t j = 0; j < 4; j++)
      uq[4 * i + j] = q[j];
  }
  return std::move(uq);
}

inline vector<double> estimate_min110(const vector<double> &uq) {
  if (0 != uq.size() % 4)
    throw std::runtime_error("'size' of input array 'uq' not "
                             "divisible by 4, 0==size%4");
  // Unit vector for 1,1,0
  const double c = 1 / std::sqrt(2.0);
  vector<double> u{c, c, 0};
  // Number of quaternions
  std::size_t n = uq.size() / 4;
  // Return value computed below
  vector<double> min110(n);
  // Cubic symmetries
  vector<double> cubic = get_symmetries();
  // Loop over random orientations
  for (std::size_t m = 0; m < n; m++) {
    std::size_t p = 4 * m;
    vector<double> r{uq[p + 0], uq[p + 1], uq[p + 2], uq[p + 3]};
    double cmax = 0.0;
    // Nested loops over symmetries
    for (std::size_t i = 0; i < 24; i++) {
      std::size_t pi = 4 * i;
      vector<double> si{cubic[pi + 0], cubic[pi + 1], cubic[pi + 2],
                        cubic[pi + 3]};
      for (std::size_t j = 0; j < 24; j++) {
        std::size_t pj = 4 * j;
        vector<double> sj{cubic[pj + 0], cubic[pj + 1], cubic[pj + 2],
                          cubic[pj + 3]};
        vector<double> sju = quaternion::rotate_vector(sj, u);
        vector<double> rsju = quaternion::rotate_vector(r, sju);
        vector<double> sirsju = quaternion::rotate_vector(si, rsju);
        double init = 0.0;
        double dot = std::fabs(
            std::inner_product(u.begin(), u.end(), sirsju.begin(), init));
        if (dot > cmax)
          cmax = dot;
      }
    }
    min110[m] = acos(cmax) * 180 / M_PI;
  }
  return std::move(min110);
}

inline vector<double> estimate_minuvw(const vector<double> &v,
                                      const vector<double> &uq) {
  if (0 != uq.size() % 4 || 3 != v.size())
    throw std::runtime_error("len(v) !=3 or 'size' of input array 'uq' not "
                             "divisible by 4, 0==size%4");
  // take a copy and normalize
  vector<double> u(v);
  {
    // normalize incoming direction v
    double nu = 1.0 / quaternion::norm(u);
    u[0] *= nu;
    u[1] *= nu;
    u[2] *= nu;
  }
  std::size_t n = uq.size() / 4;
  vector<double> minuvw(n);
  vector<double> cubic = get_symmetries();
  for (std::size_t i = 0; i < n; i++) {
    std::size_t p = 4 * i;
    vector<double> q{uq[p + 0], uq[p + 1], uq[p + 2], uq[p + 3]};
    double cmax = 0.0;
    for (std::size_t m = 0; m < 24; m++) {
      std::size_t p = 4 * m;
      vector<double> f{cubic[p + 0], cubic[p + 1], cubic[p + 2], cubic[p + 3]};
      vector<double> qs = quaternion::composition(q, f);
      quaternion::check_validity(qs);
      vector<double> r = quaternion::rotate_vector(qs, u);
      // take absolute value of each element in r
      std::transform(r.begin(), r.end(), r.begin(),
                     [](double &d) { return std::fabs(d); });
      // sort r ascending
      std::sort(r.begin(), r.end());
      double init = 0.0;
      double c = std::inner_product(u.begin(), u.end(), r.begin(), init);
      if (c > cmax)
        cmax = c;
    }
    minuvw[i] = acos(cmax) * 180 / M_PI;
  }
  return std::move(minuvw);
}

inline vector<double> estimate_min100(const vector<double> &uq) {
  /**
   * Input set of sample unit quaternions 'uq'.
   * Returns min100 distribution matching MacKenzie published result.
   */
  if (0 != uq.size() % 4)
    throw std::runtime_error(
        "'size' input array 'uq' must be divisible for 4, 0==size%4");
  std::size_t n = uq.size() / 4;
  vector<double> min100(n);
  vector<double> cubic = get_symmetries();

  // Ensure each quaternion is a unit quaternion
  for (std::size_t i = 0; i < n; i++) {
    double max = 0.0;
    std::size_t p = 4 * i;
    vector<double> r{uq[p + 0], uq[p + 1], uq[p + 2], uq[p + 3]};
    quaternion::check_validity(r);
    for (std::size_t j = 0; j < 24; j++) {
      std::size_t w = 4 * j;
      vector<double> s{cubic[w + 0], cubic[w + 1], cubic[w + 2], cubic[w + 3]};
      vector<double> t = quaternion::composition(s, r);
      vector<double> m = quaternion::to_rotation_matrix(t);
      double m1 = *std::max_element(m.begin(), m.end());
      if (m1 > max)
        max = m1;
      min100[i] = acos(max) * 180 / M_PI;
    }
  }
  return std::move(min100);
}

} // namespace CUBIC

} // namespace SPPARKS_NS

#endif
