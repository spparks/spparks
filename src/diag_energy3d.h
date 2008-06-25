/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_ENERGY3D_H
#define DIAG_ENERGY3D_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagEnergy3d : public Diag {

 public:
  DiagEnergy3d(class SPPARKS *, int, char **);
  virtual ~DiagEnergy3d();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  class AppLattice3d *applattice3d;

  int nx_local,ny_local,nz_local;

  double energy;

  FILE* fp;

};

}

#endif
