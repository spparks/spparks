/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
                NTESS, the U.S. Government retains certain rights in this software.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts/weld/jom,AppPottsWeldJOM)

#else

#ifndef SPK_APP_POTTSWELDJOM_H
#define SPK_APP_POTTSWELDJOM_H

#include "app_potts.h"
#include <stdlib.h>

namespace SPPARKS_NS {

class AppPottsWeldJOM : public AppPotts {
 public:
  AppPottsWeldJOM(class SPPARKS *, int, char **);
  virtual ~AppPottsWeldJOM() {}
  virtual void grow_app();
  virtual void init_app();
  virtual void input_app(char *, int, char **);
  virtual void site_event_rejection(int, class RandomPark *);
	
 protected:
  int nx;
  int ny; 
  int nz;
  double *MobilityOut;
  int Wwidth;
  int Wlength;
  int Wcap;
  int Haz;
  int StartWeld;
  double vel;
  int deep_width;
  int deep_length;
  int ellipsoid_depth;
  double exp_factor;
  int weld_type;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
