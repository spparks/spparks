/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "const_node.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

ConstNode::ConstNode() : Node()
{
  type = CONSTANT;
}

/* ---------------------------------------------------------------------- */

ConstNode::~ConstNode()
{
}

/* ---------------------------------------------------------------------- */

void ConstNode::write(FILE * dest)
{
  fprintf(dest,"%g",value);
}

/* ---------------------------------------------------------------------- */

void ConstNode::write_tex(FILE * dest)
{
  fprintf(dest,"%g",value);
}

/* ---------------------------------------------------------------------- */

void ConstNode::write_stack(FILE * dest)
{
  fprintf(dest," %g ",value);
}
