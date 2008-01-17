/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef NODE_H
#define NODE_H

#include "random_park.h"
#include <cmath>

namespace SPPARKS {


#define PLUS         0
#define MINUS        1
#define TIMES        2
#define DIVIDE       3
#define POW          4
#define CONSTANT     5
#define VARIABLE     6
#define VAR_DBL     7
#define VAR_INT     8
#define VAR_BOOL    9

  class Node {
  public:
    
    Node();
    virtual ~Node(){}
    
    // pure virtual functions, must be defined in child class
    
    virtual void write(FILE *) = 0;
    virtual void write_tex(FILE *) = 0;
    virtual void write_stack(FILE *) = 0;
    virtual void buffer(char **, int&) = 0;
    
    virtual double go(double *) = 0;
    virtual double go(int);
    virtual void set_left_child(Node *in) = 0;
    virtual void set_right_child(Node *in)= 0;
    virtual void clear() = 0;
    virtual void clear_ch() = 0;
    
    virtual bool equals(Node *in) = 0;
    
    
    int type;
    
};
  
}

#endif
