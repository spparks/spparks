#include "expression.h"
#include <iostream>
#include "random_park.h"
#include "error.h"
#include "state.h"
#include "tree.h"
#include "node.h"

using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

Xpression::Xpression()
{
  tree = NULL;
}

/* ---------------------------------------------------------------------- */

Xpression::~Xpression()
{

}
/* ---------------------------------------------------------------------- */
char *Xpression::get_name()
{ 
  return name;
}
/* ---------------------------------------------------------------------- */
void Xpression::set_name(char *tname)
{
  int n = strlen(tname) + 1;
  char name[n];
  strcpy(name,tname);
}
/* ---------------------------------------------------------------------- */
void Xpression::set_expression(int narg, char **arg)
{
   set_name(arg[0]);

   char inpln[80];
   strcpy(inpln,arg[1]);
   for(int t = 2; t< narg; t++) strcat(inpln,arg[t]);

   //   cout << "resulting line "<<inpln<<endl;
   if(tree)
     root = tree->from_string(inpln);
   else
     cout << "tree not set for the expression"<<endl;

}
/* ---------------------------------------------------------------------- */
double Xpression::eval(int indx)
{
  return root->go(indx);
}
