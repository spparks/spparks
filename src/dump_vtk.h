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
