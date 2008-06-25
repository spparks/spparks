/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DBL_VAR_NODE_H
#define DBL_VAR_NODE_H

#include "var_node.h"

namespace SPPARKS {

  class DblVarNode : public VarNode {
  public:
    DblVarNode();
    ~DblVarNode(){};
    inline double go(int index_in){return data_p[index_in];}
    inline double go(int index_in, int index_nb){
      //      std::cout << "nb_flag = "<<nb_flag<<std::endl;
      return (nb_flag == 0)?(data_p[index_in]):(data_p[index_nb]);
    }
    void set_data_pointer(void *in);
    
    double *data_p;
  };

}

#endif
