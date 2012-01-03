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

#ifndef SPK_LATTICE_H
#define SPK_LATTICE_H

#include "pointers.h"

namespace SPPARKS_NS {

class Lattice : protected Pointers {
 public:
  int style;                           // enum list of NONE,SC,FCC,etc
  double xlattice,ylattice,zlattice;   // lattice scale factors in 3 dims
  double a1[3],a2[3],a3[3];            // edge vectors of unit cell
  int nbasis;                          // # of basis atoms in unit cell
  double **basis;                      // fractional coords of each basis atom
                                       // within unit cell (0 <= coord < 1)
  tagint nrandom;                      // # of sites for random lattices
  double cutoff;                       // neighbor cutoff for random lattices

  Lattice(class SPPARKS *, int, char **);
  ~Lattice();
  int ncolors(int);
  int id2color(tagint, int);

private:
  double latconst;                     // lattice constant
  double origin[3];                    // lattice origin
  int orientx[3];                      // lattice orientation vecs
  int orienty[3];                      // orientx = what lattice dir lies
  int orientz[3];                      //           along x dim in box

  void add_basis(double, double, double);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Lattice style does not match dimension

Self-explanatory.

E: Cannot use coloring without domain nx,ny,nz defined

UNDOCUMENTED

E: Color stencil is incommensurate with lattice size

Since coloring induces a pattern of colors, this pattern
must fit an integer number of times into a periodic lattice.

*/
