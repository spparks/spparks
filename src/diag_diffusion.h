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
DiagStyle(diffusion,DiagDiffusion)

#else

#ifndef SPK_DIAG_DIFFUSION_H
#define SPK_DIAG_DIFFUSION_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagDiffusion : public Diag {

 public:
  DiagDiffusion(class SPPARKS *, int, char **);
  ~DiagDiffusion() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppDiffusion *appdiff;
  double deposit_success,deposit_failed;
  double nfirst_all,nsecond_all;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag_style diffusion requires app_style diffusion

Self-explanatory.

*/
