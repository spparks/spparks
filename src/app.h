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

#ifndef APP_H
#define APP_H

#include "pointers.h"

namespace SPPARKS_NS {

class App : protected Pointers {
 public:
  enum AppClasses {GENERAL,LATTICE};
  int appclass;           // one of the enum values
  char *style;            // style name of app
  double time;            // elapsed time due to events
  double time_eps;        // epsilon to avoid time precision issues

  App(class SPPARKS *, int, char **);
  virtual ~App();
  
  // pure virtual functions, must be defined in child class
  
  virtual void input(char *, int, char **) = 0;
  virtual void init() = 0;
  virtual void run(int, char **) = 0;

  // virtual functions, may be overridden in child class

  virtual void stats(char *strtmp) {strtmp[0] = '\0';};
  virtual void stats_header(char *strtmp) {strtmp[0] = '\0';};
  virtual void *extract(char *) {return NULL;};
};

}

#endif
