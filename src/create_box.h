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

#ifdef COMMAND_CLASS
CommandStyle(create_box,CreateBox)

#else

#ifndef SPK_CREATE_BOX_H
#define SPK_CREATE_BOX_H

#include "pointers.h"

namespace SPPARKS_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class SPPARKS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_box command before app_style set

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot create box with this application style

This application does not support spatial domains.

E: Cannot create box after simulation box is defined

Self-explanatory.

E: Cannot run 2d simulation with nonperiodic Z dimension

UNDOCUMENTED

E: Cannot run 1d simulation with nonperiodic Y or Z dimension

UNDOCUMENTED

E: Create_box region ID does not exist

Self-explanatory.

E: Create_box region must be of type inside

Self-explanatory.

*/
