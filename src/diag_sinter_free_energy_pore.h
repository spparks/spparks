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

#ifdef DIAG_CLASS
DiagStyle(sinter_free_energy_pore,DiagSinterFreeEnergyPore)

#else

#ifndef SPK_DIAG_SINTER_FREE_ENERGY_PORE_H
#define SPK_DIAG_SINTER_FREE_ENERGY_PORE_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagSinterFreeEnergyPore : public Diag {

 public:
  DiagSinterFreeEnergyPore(class SPPARKS *, int, char **);
  ~DiagSinterFreeEnergyPore() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);
	
protected:
  void initialize_parameters_calculation();

 private:
  class AppSinter *appsinter;
  int nlocal;
  double interfacialFE;
  bool init_flag;
	
  int xstart_, xend_;
  int ystart_, yend_;
  int zstart_, zend_;
	
	
};

}

#endif

#endif