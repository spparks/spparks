/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(sinter_free_energy,DiagSinterFreeEnergy)

#else

#ifndef SPK_DIAG_SINTER_FREE_ENERGY_H
#define SPK_DIAG_SINTER_FREE_ENERGY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagSinterFreeEnergy : public Diag {

 public:
  DiagSinterFreeEnergy(class SPPARKS *, int, char **);
  ~DiagSinterFreeEnergy() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppSinter *appsinter;
  int nlocal;
//  double density;
  double interfacialFE;
};

}

#endif

#endif
