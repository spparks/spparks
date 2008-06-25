/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "int_var_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

IntVarNode::IntVarNode() : VarNode()
{
  type = VAR_INT;
}

/* ---------------------------------------------------------------------- */

// double IntVarNode::go(int index_in)
// {
//   //  cout <<"int: "<<data_p[index_in]<<endl;
//   return data_p[index_in];
// }

/* ---------------------------------------------------------------------- */

void IntVarNode::set_data_pointer(void *in)
{
  data_p = (int *)(in);
}
