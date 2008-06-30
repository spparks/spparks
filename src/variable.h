/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef VARIABLE_H
#define VARIABLE_H

#include "sysptr.h"

namespace SPPARKS_NS {

class Variable : protected SysPtr {
 public:
  explicit Variable(class SPPARKS *);
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
