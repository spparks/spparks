/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "divide_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DivideNode::DivideNode() : Node()
{ 
  type = DIVIDE;
  left_child = NULL;
  right_child = NULL;
}

/* ---------------------------------------------------------------------- */

DivideNode::~DivideNode()
{
//   delete left_child;
//   delete right_child;
}

/* ---------------------------------------------------------------------- */

void DivideNode::write_tex(FILE *dest)
{
  fprintf(dest,"\\left( \\frac{");
  left_child->write_tex(dest);
  fprintf(dest,"}{");
  right_child->write_tex(dest);
  fprintf(dest,"}\\right)");
}

/* ---------------------------------------------------------------------- */

void DivideNode::write(FILE *dest)
{
  fprintf(dest,"(");
  left_child->write(dest);
  fprintf(dest,"/");
  right_child->write(dest);
  fprintf(dest,")");
}

/* ---------------------------------------------------------------------- */

void DivideNode::write_stack(FILE *dest)
{
  fprintf(dest,"/");
}

/* ---------------------------------------------------------------------- */

double DivideNode::go(double * var)
{
  return left_child->go(var)/right_child->go(var);
}
