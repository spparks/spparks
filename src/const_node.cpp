
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "const_node.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

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
/* ---------------------------------------------------------------------- */
void ConstNode::write(FILE * dest)
{
  fprintf(dest,"%g",value);
}

/* ---------------------------------------------------------------------- */
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
