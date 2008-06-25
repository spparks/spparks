
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "times_node.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;
/* ---------------------------------------------------------------------- */
TimesNode::TimesNode() : Node()
{ 
  type = TIMES;
  left_child = NULL;
  right_child = NULL;
}
/* ---------------------------------------------------------------------- */
TimesNode::~TimesNode()
{
//   delete left_child;
//   delete right_child;
}
/* ---------------------------------------------------------------------- */
void TimesNode::write_tex(FILE * dest)
{
  fprintf(dest,"\\left(");
  left_child->write_tex(dest);
  fprintf(dest,"*");
  right_child->write_tex(dest);
  fprintf(dest,"\\right)");
}
/* ---------------------------------------------------------------------- */
void TimesNode::write(FILE * dest)
{
  fprintf(dest,"(");
  left_child->write(dest);
  fprintf(dest,"*");
  right_child->write(dest);
  fprintf(dest,")");
}
/* ---------------------------------------------------------------------- */
void TimesNode::write_stack(FILE * dest)
{
  fprintf(dest,"*");
}
/* ---------------------------------------------------------------------- */
double TimesNode::go(double * var)
{
  return left_child->go(var) * right_child->go(var);
}
/* ---------------------------------------------------------------------- */
