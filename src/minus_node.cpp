
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "minus_node.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"
#include <iostream>

using namespace std;
using namespace SPPARKS;
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
/* ---------------------------------------------------------------------- */
