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

#ifdef DUMP_CLASS

DumpStyle(vtk,DumpVTK)

#else

#ifndef SPK_DUMP_VTK_H
#define SPK_DUMP_VTK_H

#include "dump_text.h"

namespace SPPARKS_NS {

class DumpVTK : public DumpText {
 public:
  DumpVTK(class SPPARKS *, int, char **);
  ~DumpVTK() {}
  void init_style();
  void write_header(bigint, double);
  void write_data(int, double *);
  void write_footer();
  int modify_param(int, char **);

 private:
  int type;
  int vtkflag;
  int nx,ny,nz,minval,maxval;     // set by dump_modify vtk
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
