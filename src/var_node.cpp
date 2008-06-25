/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "var_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

VarNode::VarNode() : Node()
{
  type = VARIABLE;
  nb_flag = 0;
}

/* ---------------------------------------------------------------------- */

VarNode::~VarNode()
{
}

/* ---------------------------------------------------------------------- */

void VarNode::write(FILE * dest)
{
  fprintf(dest,"var[%d]",var_index);
  // fprintf(dest,"2.5");
}

/* ---------------------------------------------------------------------- */

void VarNode::write_tex(FILE * dest)
{
  fprintf(dest,"var[%d]",var_index);
}

/* ---------------------------------------------------------------------- */

void VarNode::write_stack(FILE * dest)
{
  fprintf(dest,"var[%d]",var_index);
}

/* ---------------------------------------------------------------------- */

double VarNode::go(double * var)
{
  return var[var_index];
}
