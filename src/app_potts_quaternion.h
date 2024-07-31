/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.

   This app author:
   John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
------------------------------------------------------------------------- */

#ifdef APP_CLASS
// clang-format off
AppStyle(potts/quaternion,AppPottsQuaternion)
// clang-format on
#else

#ifndef SPK_APP_POTTS_QUATERNION_H
#define SPK_APP_POTTS_QUATERNION_H

#include "app_potts.h"
#include <array>
#include <vector>
using std::array;
using std::vector;

namespace SPPARKS_NS {

class AppPottsQuaternion : public AppPotts {
public:
  AppPottsQuaternion(class SPPARKS *, int, char **);
  virtual ~AppPottsQuaternion() {}
  virtual void init_app();
  virtual void grow_app();
  virtual double site_energy(int i);
  virtual void site_event_rejection(int, class RandomPark *);
  virtual double site_propensity(int);

  friend struct SiteState;
  struct SiteState {
  public:
    SiteState(int s, const array<double, 4> &qs) : spin(s), q(qs) {}
    int spin;
    array<double, 4> q;
  };

private:
  // symmetries: expect 'cubic'(24) or 'hcp'(12)
  vector<double> symmetries;
  // Possible extension to multiple symmetry types;
  // Each site carries a symmetry type; then symmetries for
  // that site can be looked up using this mapping.
  // map<key=string,vector<double>> symmetries;

  // orientation stored as quaternion
  double *q0, *qx, *qy, *qz;
  int *unique_neigh;
  void flip_site(int site, const SiteState &s);
};

} // namespace SPPARKS_NS

#endif
#endif

    /* ERROR/WARNING messages:

    E: Illegal ... command

    Self-explanatory.  Check the input script syntax and compare to the
    documentation for the command.  You can use -echo screen as a
    command-line option when running SPPARKS to see the offending
    line.

    E: One or more sites have invalid values

    The application only allows sites to be initialized with specific
    values.

    */
