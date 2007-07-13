
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "divide_node.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"
#include <iostream>

using namespace std;
using namespace SPPARKS;
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
/* ---------------------------------------------------------------------- */
