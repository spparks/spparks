/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef EXPRESSION_H
#define EXPRESSION_H

#include "tree.h"

namespace SPPARKS {
  
class Xpression {
    
  protected:
    char name[50];
    class Tree *tree;
    class Node *root;
    
  public:
    Xpression();
    ~Xpression();
    
    char *get_name();
    void set_name(char *);
    void set_expression(int, char**);

    inline void set_tree(class Tree *t){tree = t;}

    double eval(int);
    double sum(int);
   
    class Neighborhood *nbr;
    int nb_flag;
    inline void set_neighborhood(class Neighborhood *nb)
      {nb_flag = 1; nbr = nb;}
  };
  
}

#endif

