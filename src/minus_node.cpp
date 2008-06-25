/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "minus_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

MinusNode::MinusNode() : Node()
{ 
  type = MINUS;
  left_child = NULL;
  right_child = NULL;
}

/* ---------------------------------------------------------------------- */

MinusNode::~MinusNode()
{
//   delete left_child;
//   delete right_child;
}

/* ---------------------------------------------------------------------- */

void MinusNode::write_tex(FILE * dest)
{
  fprintf(dest,"\\left(");
  left_child->write_tex(dest);
  fprintf(dest,"-");
  right_child->write_tex(dest);
  fprintf(dest,"\\right)");
}

/* ---------------------------------------------------------------------- */

void MinusNode::write(FILE * dest)
{
  fprintf(dest,"(");
  left_child->write(dest);
  fprintf(dest,"-");
  right_child->write(dest);
  fprintf(dest,")");
}

/* ---------------------------------------------------------------------- */

void MinusNode::write_stack(FILE * dest)
{
  fprintf(dest,"-");
}

/* ---------------------------------------------------------------------- */

double MinusNode::go(double *var)
{
  return left_child->go(var) - right_child->go(var);
}
