/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "pow_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

PowNode::PowNode() : Node()
{ 
  type = POW;
  left_child = NULL;
  right_child = NULL;
}

/* ---------------------------------------------------------------------- */

PowNode::~PowNode()
{
//   delete left_child;
//   delete right_child;
}

/* ---------------------------------------------------------------------- */

void PowNode::write_tex(FILE *dest)
{
  fprintf(dest,"\\left( {");
  left_child->write_tex(dest);
  fprintf(dest,"} \\right) ^{ \\left(");
  right_child->write_tex(dest);
  fprintf(dest,"\\right) }");
}

/* ---------------------------------------------------------------------- */

void PowNode::write(FILE *dest)
{
  fprintf(dest,"(");
  left_child->write(dest);
  fprintf(dest,")**(");
  right_child->write(dest);
  fprintf(dest,")");
}

/* ---------------------------------------------------------------------- */

void PowNode::write_stack(FILE *dest)
{
  fprintf(dest,"^");
}

/* ---------------------------------------------------------------------- */

double PowNode::go(double * var)
{
  return pow(left_child->go(var),right_child->go(var));
}
