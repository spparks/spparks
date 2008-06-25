/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef CONST_NODE_H
#define CONST_NODE_H

#include "random_park.h"
#include <cmath>
#include "node.h"
#include <iostream>

using namespace std;

namespace SPPARKS {

class ConstNode : public Node {
 public:
  ConstNode();
  ~ConstNode();


  void write(FILE *);
  void write_tex(FILE *);
  void write_stack(FILE *);
  inline void buffer(char **buf, int &n){
    sprintf(buf[n],"%d %lg ",CONSTANT, value); n++;
  };

  inline double go(double * var){return value;}
  inline void set_value(double value_in){value = value_in;}
  inline void set_left_child(Node *in){};
  inline void set_right_child(Node *in){};
  inline void clear(){delete this;};
  inline void clear_ch(){};
  inline double get_value() {return value;}

  private:

  double value;

};

}

#endif
