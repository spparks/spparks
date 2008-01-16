/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef DBL_VAR_NODE_H
#define DBL_VAR_NODE_H

#include <cmath>
#include "node.h"
#include <iostream>
#include "var_node.h"

namespace SPPARKS {

  class DblVarNode : public VarNode {
  public:
    DblVarNode();
    ~DblVarNode(){};
    inline double go(int index_in){return data_p[index_in];}
    
    void set_data_pointer(void *in);
    
    double *data_p;
  };
}


#endif