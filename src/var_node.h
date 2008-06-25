/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef VAR_NODE_H
#define VAR_NODE_H

#include "stdio.h"
#include "string.h"
#include "node.h"

namespace SPPARKS_NS {

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

  inline bool equals(Node *in){
    if (in->type == type)
      if(static_cast<VarNode*>(in)->get_value()==var_index)
	return true;

    return false;
  };

  char var_name[16];
  inline void set_name(char *name_in){strcpy(var_name,name_in);}
  inline char *get_name(){return var_name;}
  inline void set_data_type(int in){data_type = in;}
  virtual void set_data_pointer(void *in){data_p = in;}
  inline void set_nb(){nb_flag = 1;}

  void *data_p;
  int nb_flag;
  private:

  int var_index;
  int data_type;
};

}

#endif
