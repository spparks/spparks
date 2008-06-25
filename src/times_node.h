/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef TIMES_NODE_H
#define TIMES_NODE_H

#include "node.h"

namespace SPPARKS_NS {

class TimesNode : public Node {
 public:
  TimesNode();
  ~TimesNode();

  double go(double *);
  inline double go(int index_in){
    return left_child->go(index_in) * right_child->go(index_in);
  }
  inline double go(int index_in, int index_nb){
    return left_child->go(index_in, index_nb) * right_child->go(index_in,index_nb);
  }
  void write(FILE *);
  void write_tex(FILE *);
  void write_stack(FILE *);
  inline void buffer(char **buf, int &n){
    sprintf(buf[n],"%d ",TIMES); n++;
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

  inline bool equals(Node *in){
    if (in->type == type){
      Node* left = static_cast<TimesNode*>(in)->get_left_child();
      Node* right = static_cast<TimesNode*>(in)->get_right_child();
      if(left_child->equals(left) && right_child->equals(right))
	return true;
    }
    return false;
  };

 private:
  Node *left_child;
  Node *right_child;
};

}

#endif
