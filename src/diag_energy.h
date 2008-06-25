/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_ENERGY_H
#define DIAG_ENERGY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagEnergy : public Diag {

 public:
  DiagEnergy(class SPPARKS *, int, char **);
  virtual ~DiagEnergy();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  class AppLattice *applattice;

  int nlocal;

  double energy;

  FILE* fp;

};

}

#endif
