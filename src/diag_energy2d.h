/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_ENERGY2D_H
#define DIAG_ENERGY2D_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS {

class DiagEnergy2d : public Diag {

 public:
  DiagEnergy2d(class SPK *, int, char **);
  virtual ~DiagEnergy2d();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  class AppLattice2d *applattice2d;

  int nx_local,ny_local;

  double energy;

  FILE* fp;

};

}

#endif
