/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef VAR_NODE_H
#define VAR_NODE_H

#include "random_park.h"
#include <cmath>
#include "node.h"
#include <iostream>

using namespace std;
namespace SPPARKS {

class VarNode : public Node {
 public:
  VarNode();
  ~VarNode();

  void write(FILE *);
  void write_tex(FILE *);
  void write_stack(FILE *);
  inline void buffer(char **buf, int &n){
    sprintf(buf[n],"%d %d ",VARIABLE, var_index); n++;
  };
  double go(double *);
  inline void set_value(int in){var_index = in;}
  inline int get_value() {return var_index;}
  inline void set_left_child(Node *in){};
  inline void set_right_child(Node *in){};
  inline void clear(){delete this;};
  inline void clear_ch(){};

  private:

  int var_index;

};

}

#endif
