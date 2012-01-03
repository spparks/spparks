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

#ifndef SPK_VARIABLE_H
#define SPK_VARIABLE_H

#include "pointers.h"

namespace SPPARKS_NS {

class Variable : protected Pointers {
 public:
  Variable(class SPPARKS *);
  ~Variable();
  void set(int, char **);
  void set(char *, char *);
  int next(int, char **);
  int find(char *);
  int equalstyle(int);
  char *retrieve(char *);
  double compute_equal(int);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables arrays can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *index;              // next available value for each variable
  char ***data;            // str value of each variable's values
  int precedence[7];       // precedence level of math operators

  double PI;

  void remove(int);
  void copy(int, char **, char **);
  double evaluate(char *);
  int find_matching_paren(char *, int, char *&);
  int math_function(char *, char *, double *, int &);
  int is_constant(char *);
  double constant(char *);
  int keyword(char *, double &);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot redefine variable as a different style

An equal-style variable can be re-defined but only if it was
originally an equal-style variable.

E: World variable count doesn't match # of partitions

A world-style variable must specify a number of values equal to the
number of processor partitions.

E: Universe/uloop variable count < # of partitions

A universe or uloop style variable must specify a number of values >= to the
number of processor partitions.

E: All universe/uloop variables must have same # of values

Self-explanatory.

E: Variable name must be alphanumeric or underscore characters

Self-explanatory.

E: Invalid variable in next command

Self-explanatory.

E: All variables in next command must be same style

Self-explanatory.

E: Invalid variable style with next command

Variable styles {equal} and {world} cannot be used in a next
command.

E: Invalid syntax in variable formula

Self-explanatory.

E: Invalid variable name in variable formula

Variable name is not recognized.

E: Invalid variable evaluation in variable formula

A variable used in a formula could not be evaluated.

E: Invalid math function in variable formula

The math function is not recognized.

E: Invalid keyword in variable formula

UNDOCUMENTED

E: Divide by 0 in variable formula

Self-explanatory.

E: Power by 0 in variable formula

Self-explanatory.

E: Sqrt of negative in variable formula

Self-explanatory.

E: Log of zero/negative in variable formula

Self-explanatory.

E: Arcsin of invalid value in variable formula

Argument of arcsin() must be between -1 and 1.

E: Arccos of invalid value in variable formula

Argument of arccos() must be between -1 and 1.

*/
