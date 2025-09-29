#ifndef BURGERSORIENATIONRELATION_BCC_TO_HCP_H
#define BURGERSORIENATIONRELATION_BCC_TO_HCP_H

#include <math.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <numbers>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "quaternion.h"

using SPPARKS_NS::quaternion::composition;
using std::asin;
using std::vector;

namespace SPPARKS_NS {

namespace BurgersOrienationRelation_BCC_to_HCP {

inline auto rot_x = [](double phi) {
  const double pi = M_PI;
  const double rad = phi * pi / 180.0 / 2.0;

  // Create a std::vector<double> with 4 elements
  std::vector<double> rot(4);

  rot[0] = std::cos(rad);
  rot[1] = std::sin(rad);
  rot[2] = 0.0;
  rot[3] = 0.0;

  return rot;
};

inline auto rot_y = [](double phi) {
  const double pi = M_PI;
  const double rad = phi * pi / 180.0 / 2.0;

  // Create a std::vector<double> with 4 elements
  std::vector<double> rot(4);

  rot[0] = std::cos(rad);
  rot[1] = 0.0;
  rot[2] = std::sin(rad);
  rot[3] = 0.0;

  return rot;
};

inline auto rot_z = [](double phi) {
  const double pi = M_PI;
  const double rad = phi * pi / 180.0 / 2.0;

  // Create a std::vector<double> with 4 elements
  std::vector<double> rot(4);

  rot[0] = std::cos(rad);
  rot[1] = 0.0;
  rot[2] = 0.0;
  rot[3] = std::sin(rad);

  return rot;
};

inline vector<double> get_variant(const vector<double> &q0, const vector<double> &q1,
                           const vector<double> &q2, const vector<double> &q3,
                           const vector<double> &q4) {
  /**

  Returns so-called active orientations.

  SPPARKS stores the active quaternion orientations on
  each lattice site.

  Active or passive orientations can be used to compute disorientations
  between variants -- provided that the disorientation calculation
  is tweaked to take into account 'active' or 'passive' variants.

  See compute_disorientation in disorientation.h

  For plotting, the output of this function
  will orient/rotate the variant and depicts
  the variant in its actual orientation.
  */

  auto r1 = composition(q0, q1);
  auto r2 = composition(r1, q2);
  auto r3 = composition(r2, q3);
  auto r4 = composition(r3, q4);
  return std::move(r4);
}

/*
 * Return 12 Variants
 * ------------------
 *
 * In titanium, alpha phase (HCP) randomly nucleates with an orientation
 * from one of these variants. Orientation of variants are relative to the
 * parent BCC orientation. The parent BCC is often called the 'prior'.
 *
 * Number of misorientations between variants:
 *     78=(12*12-12)/2 + 12 = 66 + 12
 *
 * Misorientation between variants should follow the expected statistics.
 */

inline vector<vector<double>>
get_variants(const vector<double> &q0 = vector<double>{1, 0, 0, 0}) {
  /**
   * Compute BOR variants.  Variants are wrt to 'prior'
   *
   * arg q0: orientation of cube
   *
   */
  const double sqrt2 = M_SQRT2;
  const double sqrt3 = std::sqrt(3.0);
  const double pi = M_PI;
  const double phi = 180 * asin(1.0 / sqrt3) / pi;
  /**
   ***************************
   * Orient hexagonal prism
   ***************************
   */
  // xy edge
  auto qxy = rot_x(0);
  // xz edge
  auto qxz = rot_x(90);
  // yz edge
  auto qyz = composition(rot_z(90), rot_x(90));

  vector<vector<double>> variants(12);
  // # Table 2 Karthikeyan (110)
  variants[0] = get_variant(q0, qxy, rot_z(-45), rot_x(90), rot_z(phi));
  variants[1] = get_variant(q0, qxy, rot_z(-45), rot_x(90), rot_z(60 - phi));

  // # Table 2 Karthikeyan (1b10) or (11b0)
  variants[2] = get_variant(q0, qxy, rot_z(45), rot_x(90), rot_z(60 - phi));
  variants[3] = get_variant(q0, qxy, rot_z(45), rot_x(90), rot_z(phi));

  // #Table 2 Karthikeyan(101)
  variants[4] = get_variant(q0, qxz, rot_z(-45), rot_x(90), rot_z(phi));
  variants[5] = get_variant(q0, qxz, rot_z(-45), rot_x(90), rot_z(60 - phi));

  // #Table 2 Karthikeyan(1b01)
  variants[6] = get_variant(q0, qxz, rot_z(45), rot_x(90), rot_z(phi));
  variants[7] = get_variant(q0, qxz, rot_z(45), rot_x(90), rot_z(60 - phi));

  // #Table 2 Karthikeyan(011)
  variants[8] = get_variant(q0, qyz, rot_z(-45), rot_x(90), rot_z(phi));
  variants[9] = get_variant(q0, qyz, rot_z(-45), rot_x(90), rot_z(60 - phi));

  // #Table 2 Karthikeyan(01b1)
  variants[10] = get_variant(q0, qyz, rot_z(45), rot_x(90), rot_z(phi));
  variants[11] = get_variant(q0, qyz, rot_z(45), rot_x(90), rot_z(60 - phi));
  return variants;
}

} // namespace BurgersOrienationRelation_BCC_to_HCP

} // namespace SPPARKS_NS

#endif // BURGERSORIENATIONRELATION_BCC_TO_HCP_H
