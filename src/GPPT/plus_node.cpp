
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "plus_node.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;
/* ---------------------------------------------------------------------- */
PlusNode::PlusNode() : Node()
{ 
  type = PLUS;
  left_child = NULL;
  right_child = NULL;
}
/* ---------------------------------------------------------------------- */
PlusNode::~PlusNode()
{
//   delete left_child;
//   delete right_child;
}
/* ---------------------------------------------------------------------- */
void PlusNode::write_tex(FILE * dest)
{
    fprintf(dest,"\\left(");
    left_child->write_tex(dest);
    fprintf(dest,"*");
    right_child->write_tex(dest);
    fprintf(dest,"\\right)");
 
}
/* ---------------------------------------------------------------------- */
void PlusNode::write(FILE * dest)
{
  fprintf(dest,"(");
  left_child->write(dest);
  fprintf(dest,"*");
  right_child->write(dest);
  fprintf(dest,")");
}
/* ---------------------------------------------------------------------- */
void PlusNode::write_stack(FILE * dest)
{
  fprintf(dest,"*");
}
/* ---------------------------------------------------------------------- */
double PlusNode::go(double * var)
{
  return left_child->go(var) + right_child->go(var);
}
/* ---------------------------------------------------------------------- */
