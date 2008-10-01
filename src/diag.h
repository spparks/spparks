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

#ifndef DIAG_H
#define DIAG_H

#include "pointers.h"

namespace SPPARKS_NS {

class Diag : protected Pointers {
 public:
  Diag(class SPPARKS *, int, char **);
  virtual ~Diag();

  int check_time(double, int);
  void setup_time(double);

  // virtual functions with default definitions
  // may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};

  // pure virtual functions, must be defined in child class
  
  virtual void init(double) = 0;
  virtual void compute(double, int, int) = 0;

 protected:
  char *style;
  int me,nprocs;

  double diag_time,diag_delta,diag_scale,diag_t0;
  int diag_nrepeat,diag_irepeat,diag_ilogfreq;
  int stats_flag;

};

}

#endif
