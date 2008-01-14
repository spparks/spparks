#ifndef XPRESSION_H
#define XPRESSION_H

#include "sysptr.h"
#include "state.h"
#include "tree.h"
#include "stdio.h"

namespace SPPARKS{
  
  class Xpression{
    
  protected:
    
    char *name;
    State *state;

    Tree *tree;
    Node *root;
    
  public:
    Xpression();
    ~Xpression();
    
    char *get_name();
    void set_name(char *);
    void set_expression(int, char**);
    inline void set_state(State * st){state = st;}
    inline void set_tree(Tree *t){tree = t;}

    double eval(int);
    
    
  };
  
}

#endif

