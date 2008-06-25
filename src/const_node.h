/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef CONST_NODE_H
#define CONST_NODE_H

#include "node.h"

namespace SPPARKS_NS {

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
  inline double go(int index_in){return value;}
  inline double go(int index_in, int index_nb){return value;}
  inline void set_value(double value_in){value = value_in;}
  inline void set_left_child(Node *in){};
  inline void set_right_child(Node *in){};
  inline void clear(){delete this;};
  inline void clear_ch(){};
  inline double get_value() {return value;}

  inline bool equals(Node *in){
    if (in->type == type)
      if(static_cast<ConstNode*>(in)->get_value()==value)
	return true;
    return false;
  };
  private:

  double value;

};

}

#endif
