#ifndef XPRESSION_H
#define XPRESSION_H

#include "sysptr.h"
#include "neighborhood.h"
#include "tree.h"
#include "stdio.h"


namespace SPPARKS{
  
  class Xpression{
    
  protected:
    
    char name[50];


    Tree *tree;
    Node *root;
    
  public:
    Xpression();
    ~Xpression();
    
    char *get_name();
    void set_name(char *);
    void set_expression(int, char**);

    inline void set_tree(Tree *t){tree = t;}

    double eval(int);
    double sum(int);
   
    Neighborhood *nbr;
    int nb_flag;
    inline void set_neighborhood(Neighborhood *nb)
      {nb_flag = 1; nbr = nb;}
    
  };
  
}

#endif

