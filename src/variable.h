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

#ifndef VARIABLE_H
#define VARIABLE_H

#include "pointers.h"

namespace SPPARKS_NS {

class Variable : protected Pointers {
 public:
  Variable(class SPPARKS *);
  ~Variable();
  void set(int, char **);
  void set(char *, char *);
  int next(int, char **);
  char *retrieve(char *);

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

  int find(char *);
  void remove(int);
  void copy(int, char **, char **);
  double evaluate(char *);
  int find_matching_paren(char *, int, char *&);
  int math_function(char *, char *, double *, int &);
};

}

#endif
