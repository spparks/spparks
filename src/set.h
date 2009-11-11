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

#ifdef COMMAND_CLASS
CommandStyle(set,Set)

#else

#include "pointers.h"

namespace SPPARKS_NS {

class Set : protected Pointers {
 public:
  Set(class SPPARKS *);
  void command(int, char **);

 private:
  int siteindex;
  int count;
  int ivalue,ivaluelo,ivaluehi;
  double dvalue,dvaluelo,dvaluehi;
  int loopflag,regionflag,iregion;
  double fraction;

  struct Condition {                     // list of if-test conditions
    int lhs,type,index,stride;
    int op;
    int irhs;
    double drhs;
  };
  Condition *cond;
  int ncondition;

  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appoff;
    
  void set_value(int);
  void set_range(int);
  void set_displace(int);
  int condition(int);
};

}

#endif
