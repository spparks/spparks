/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "string.h"
#include "expression.h"
#include "node.h"
#include "neighborhood.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Expression::Expression()
{
  tree = NULL;
  memset(name,0,50);
}

/* ---------------------------------------------------------------------- */

Expression::~Expression()
{
}

/* ---------------------------------------------------------------------- */

char *Expression::get_name()
{ 
  return &name[0];
}

/* ---------------------------------------------------------------------- */

void Expression::set_name(char *tname)
{
  strcpy(name,tname);
}

/* ---------------------------------------------------------------------- */

void Expression::set_expression(int narg, char **arg)
{
   set_name(arg[0]);

   char inpln[80];
   strcpy(inpln,arg[1]);
   for(int t = 2; t< narg; t++) strcat(inpln,arg[t]);
  
   //  cout << "resulting line "<<inpln<<endl;
   root = tree->from_string(inpln);
}

/* ---------------------------------------------------------------------- */

double Expression::eval(int indx)
{
  return root->go(indx);
}

/* ---------------------------------------------------------------------- */

double Expression::sum(int ndx)
{
  double sum = 0.0;
  int ndx_nb;
//   cout<<"sum invoked with index "<<ndx<<endl;
//   cout <<"neighbor number at this index is ";
//   cout <<(int)nbr->numneigh[ndx]<<endl;
  int num_nbr = nbr->numneigh[ndx];
  int *nbrs = nbr->neighbor[ndx];
  for(int n = 0; n < num_nbr; n++){
    sum += root->go(ndx, nbrs[n]);
  }
  return sum;
}
