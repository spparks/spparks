/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef POW_NODE_H
#define POW_NODE_H

#include "random_park.h"
#include <cmath>
#include "node.h"
#include <iostream>

using namespace std;
namespace SPPARKS {

class PowNode : public Node {
 public:
  PowNode();
  ~PowNode();

  double go(double *);

  void write(FILE *);
  void write_tex(FILE *);
  void write_stack(FILE *);
  inline void buffer(char **buf, int &n){
    sprintf(buf[n],"%d ",POW); n++;
    left_child->buffer(buf, n);
    right_child->buffer(buf, n);
  };
  inline Node* get_left_child(){return left_child;};
  inline Node* get_right_child(){return right_child;};
  inline void set_left_child(Node *in){left_child = in;};
  inline void set_right_child(Node *in){right_child = in;};
  inline void clear()
    {
      left_child->clear(); 
      right_child->clear();
      delete this;
    }
  inline void clear_ch()
    {
      left_child->clear(); 
      right_child->clear();
    }
  private:

  Node *left_child;
  Node *right_child;

};

}

#endif
