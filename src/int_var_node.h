/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef INT_VAR_NODE_H
#define INT_VAR_NODE_H

#include "var_node.h"

namespace SPPARKS_NS {

class IntVarNode : public VarNode {
  public:
    IntVarNode();
    ~IntVarNode(){};
    inline double go(int index_in){return data_p[index_in];}
    inline double go(int index_in, int index_nb){
      return (nb_flag == 0)?(data_p[index_in]):(data_p[index_nb]);
    }
    
    void set_data_pointer(void *in);
    
    int *data_p;
  };

}

#endif
