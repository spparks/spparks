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

#ifdef DIAG_CLASS
DiagStyle(sinter_density,DiagSinterDensity)

#else

#ifndef SPK_DIAG_SINTER_DENSITY_H
#define SPK_DIAG_SINTER_DENSITY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagSinterDensity : public Diag {

 public:
  DiagSinterDensity(class SPPARKS *, int, char **);
  ~DiagSinterDensity() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);
  
 protected:
  void initialize_parameters_density_calculation(); 

 private:
  class AppSinter *appsinter;
  int nlocal;
  double density;
  
  int xstart_density, xend_density;
  int ystart_density, yend_density;
  int zstart_density, zend_density;
};

}

#endif

#endif
