/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef NODE_H
#define NODE_H

#include "random_park.h"
#include <cmath>
#include "sysptr.h"

namespace SPPARKS {


#define PLUS         0
#define MINUS        1
#define TIMES        2
#define DIVIDE       3
#define POW          4
#define CONSTANT     5
#define VARIABLE     6

class Node : protected SysPtr {
 public:

  Node();
  virtual ~Node(){}

  // pure virtual functions, must be defined in child class

  virtual void write(FILE *) = 0;
  virtual void write_tex(FILE *) = 0;
  virtual void write_stack(FILE *) = 0;
  virtual void buffer(char **, int&) = 0;

  virtual double go(double *) = 0;
  virtual void set_left_child(Node *in) = 0;
  virtual void set_right_child(Node *in)= 0;
  virtual void clear() = 0;
  virtual void clear_ch() = 0;


  int type;

};

}

#endif
