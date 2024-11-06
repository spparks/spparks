/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.

   hcp_symmetries.h
   (Note: the acronym 'hcp' stands for 'hexagonal close-packed')
   Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
   Author: Meg McCarthy, megmcca@sandia.gov, Sandia National Laboratories
   Date: March 8 2024

------------------------------------------------------------------------- */
#ifndef SPK_HCP_SYMMETRIES_H
#define SPK_HCP_SYMMETRIES_H

#include "quaternion.h"
#include <algorithm>
#include <cstdio>
#include <math.h>
#include <vector>
using std::vector;

namespace SPPARKS_NS {

namespace HCP {

inline vector<double> get_symmetries() {
  // see also https://doi.org/10.1107/S0108767391006864
  // Hexagonal rotational symmetries defined using 'axis-angle' unit
  // quaternion representation.
  // Rotation angle theta, axis n.
  // q=(cos(theta/2),sin(theta/2)*[nx,ny,nz])
  //    where n = nx,ny,nz is a unit vector and is referred to as 'spindle'.
  // Note: Rodrigues vector (3-component version of quaternion):
  // d := (1/q_0)(q_1; q_2; q_3) = tan(theta/2)*n
  //    where tan(theta/2) = sin(theta/2)/cos(theta/2)
  //    and n is the same spindle defined above.
  //
  // # Constants used to define unit quaternions
  // # urX: unit radian for degrees X
  // # ucX: output cosine for degrees X
  // # usX: output sine for degrees X
  const double b = 1.0 / 2;
  const double c = 1 / std::sqrt(2.0);
  const double pi = atan(1) * 4; // helper
  const double ur60 = pi / 3;
  const double ur120 = 2.0 * pi / 3;
  const double ur180 = pi;

  const double uc30 = std::sqrt(3.0) / 2;
  const double us30 = 1.0 / 2;
  const double uc60 = 1.0 / 2;
  const double us60 = std::sqrt(3.0) / 2;
  const double uc90 = 0;
  const double us90 = 1.0;
  const double uc120 = -1.0 / 2;
  const double us120 = std::sqrt(3.0) / 2;
  const double uc150 = -std::sqrt(3.0) / 2;
  const double us150 = 1.0 / 2;
  const double uc180 = -1.0;
  const double us180 = 0;

  // # Identity (no rotation or equivalent full rotation for axis)
  const vector<double> q1{1, 0, 0, 0};

  // #####################################################
  // ### Conventional cell HCP has two 'spindle' locations:
  // ### 1: Parallel to unit axis c (long axis), centered
  // ###    - One 6-fold axis of rotation (60dg/rotation = 720dg total)
  // ###    - Through "top/bottom" of HCP unit cell
  // ###    - Adds 6 symmetry operations (1 axis x 6 rotations/axis)
  // ### 2: Perpendicular to unit axis c (long axis), in horizontal mirror plane
  // ###    - Three 2-fold axes of rotation (180dg/rotation = 360dg total)
  // ###    - Cut "faces" (vertical mirror 1) and "corners" (vertical mirror 2)
  // of HCP unit cell
  // ###    - Adds 6 symmetry operations (2 v. mirror planes x 3
  // orientations/plane x (2 rotations/axis - 1 symmetric operation in
  // identity))
  // ##

  // ##
  // # z-axis spindle (z || c)
  // # +-60 degrees, +-120 degrees, 180 degrees
  const vector<double> q2{uc30, 0, 0, us30};
  const vector<double> q3{uc60, 0, 0, us60};
  const vector<double> q4{uc30, 0, 0, -us30};
  const vector<double> q5{uc60, 0, 0, -us60};
  const vector<double> q6{uc90, 0, 0, us90};

  // ##
  // # 'corner' spindles (x || a_1)
  // # all are 180 degrees
  const vector<double> q7{uc90, us90, 0, 0};
  const vector<double> q8{uc90, uc60 * us90, us60 * us90, 0};
  const vector<double> q9{uc90, uc120 * us90, us120 * us90, 0};

  // ##
  // # 'face' spindles (x,y @ 30 degrees)
  // # all are 180 degrees
  const vector<double> q10{uc90, uc30 * us90, us30 * us90, 0};
  const vector<double> q11{uc90, uc90 * us90, us90 * us90, 0};
  const vector<double> q12{uc90, uc150 * us90, us150 * us90, 0};

  vector<vector<double>> symmetries(12);
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

  // Copy symmetries to contiguous vector
  vector<double> uq(12 * 4);
  for (std::size_t i = 0; i < 12; i++) {
    vector<double> q = symmetries[i];
    for (std::size_t j = 0; j < 4; j++)
      uq[4 * i + j] = q[j];
  }
  return std::move(uq);
}

} // namespace HCP

} // namespace SPPARKS_NS

#endif
