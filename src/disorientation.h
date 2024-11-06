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
#ifndef SPK_DISORIENTATION_H
#define SPK_DISORIENTATION_H

#include "quaternion.h"
#include <algorithm>
#include <cstdio>
#include <math.h>
#include <vector>
using std::vector;

namespace SPPARKS_NS {

namespace disorientation {

inline double compute_disorientation(const vector<double> &symmetries,
                                     const vector<double> &q1,
                                     const vector<double> &q2) {
  // Return disorientation scalar angle between input quaternions q1 and q2
  //  on the the basis of input symmetries;
  if (0 != symmetries.size() % 4)
    throw std::runtime_error(
        "'size' input array 'symmetries' must be divisible for 4, 0==size%4");
  double t = 0.0;
  std::size_t num_symmetries = symmetries.size() / 4;
  for (std::size_t i = 0; i < num_symmetries; i++) {
    std::size_t p = 4 * i;
    vector<double> qr{symmetries[p + 0], symmetries[p + 1], symmetries[p + 2],
                      symmetries[p + 3]};
    vector<double> q = quaternion::composition(qr, q1);
    vector<double> qc = quaternion::conjugate(q);
    vector<double> m = quaternion::composition(q2, qc);
    double trace = quaternion::trace_of_rotation_matrix(m);
    if (trace > t)
      t = trace;
  }
  double di = acos(0.5 * (t - 1)) * 180 / M_PI;
  return di;
}

inline vector<double>
estimate_disorientation_distribution(const vector<double> &uq,
                                     const vector<double> &symmetries) {
  /**
   * Input set of sample unit quaternions 'uq'.
   * Returns distribution of disorientations matching MacKenzie published
   * distribution.
   */
  if (0 != uq.size() % 4)
    throw std::runtime_error(
        "'size' input array 'uq' must be divisible for 4, 0==size%4");
  if (0 != symmetries.size() % 4)
    throw std::runtime_error(
        "'size' input array 'symmetries' must be divisible for 4, 0==size%4");
  // Reference aka identity quaternion
  vector<double> qref{1.0, 0.0, 0.0, 0.0};

  std::size_t n = uq.size() / 4;
  vector<double> di(n);

  // Ensure each quaternion is a unit quaternion
  for (std::size_t i = 0; i < n; i++) {
    std::size_t p = 4 * i;
    vector<double> q{uq[p + 0], uq[p + 1], uq[p + 2], uq[p + 3]};
    quaternion::check_validity(q);
    di[i] = compute_disorientation(symmetries, qref, q);
  }
  return std::move(di);
}

} // namespace disorientation

} // namespace SPPARKS_NS

#endif
