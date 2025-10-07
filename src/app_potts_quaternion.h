/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
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
  // Possible extension to handle multiple phases/symmetry types.
  // Each site carries a key for the symmetry type: integer;
  // symmetries for that site can be looked up using this mapping.
  // map<key=integer,vector<double>> symmetries;

  // site orientation stored as quaternion
  double *q0, *qx, *qy, *qz;
  // Read-Shockley disorientation cutoff
  double theta_cut;
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
command-line option when running SPPARKS to see the offending line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
