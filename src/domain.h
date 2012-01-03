/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_DOMAIN_H
#define SPK_DOMAIN_H

#include "pointers.h"

namespace SPPARKS_NS {

class Domain : protected Pointers {
 public:
  int me,nprocs;                    // proc info
  int procgrid[3];                  // assigned # of procs in each dim
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs

  int box_exist;                    // 0 = not yet created, 1 = exists
  int dimension;                    // 1,2,3
  int nonperiodic;                  // 0 = periodic in all dims
                                    // 1 = non-periodic in any dim
  int xperiodic,yperiodic,zperiodic;  // 0 = non-periodic, 1 = periodic
  int periodicity[3];               // xyz periodicity as array

  double xprd,yprd,zprd;                               // global domain
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // global box bounds
  double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;    // my portion of box
  int nx,ny,nz;               // size of box in lattice spacings
                              // only set for periodic sites from create_sites

  class Lattice *lattice;                  // user-defined lattice

  int nregion;                             // # of defined Regions
  int maxregion;                           // max # list can hold
  class Region **regions;                  // list of defined Regions

  Domain(class SPPARKS *);
  ~Domain();

  void set_box();
  void set_lattice(int, char **);
  void add_region(int, char **);
  int find_region(char *);
  void set_boundary(int, char **);

  void procs2domain_1d();
  void procs2domain_2d();
  void procs2domain_3d();
};

}

#endif

/* ERROR/WARNING messages:

E: Box bounds are invalid

Lo bound >= hi bound.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Reuse of region ID

Self-explanatory.

E: Invalid region style

Self-explanatory.

E: Boundary command currently only supported by on-lattice apps

UNDOCUMENTED

E: App style proc count is not valid for 1d simulation

There can only be 1 proc in y and z dimensions for 1d models.

E: App style proc count is not valid for 2d simulation

There can only be 1 proc in the z dimension for 2d models.

*/
