/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef VARIABLE_H
#define VARIABLE_H

#include "sysptr.h"

namespace SPPARKS {

class Variable : protected SysPtr {
 public:
  explicit Variable(class SPK *);
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

  int find(char *);
  void copy(int, char **, char **);
  char *evaluate(char *);
  void remove(int);

 private:
  Variable(); // Not a sane operation.
  Variable(const Variable&); // Not a sane operation.
  Variable& operator=(const Variable&); // Not a sane operation.
};

}

#endif
