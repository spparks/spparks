/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef TREE_H
#define TREE_H

#include "node.h"
#include <cmath>
#include "random_park.h"
#include "state.h"
#include "sysptr.h"


namespace SPPARKS {
  
#define MUT_CONST    0
#define MUT_VAR      1
#define MUT_OPER     2
#define MUT_BRANCH   3
#define MUT_CROSS    4
#define MUT_SWAP     5
  
  class Tree{
  public:
    
    Tree();
    ~Tree();
    
    Node *root;
    
    Node *build_tree(RandomPark*, int);
    
    Node *random_node(RandomPark*, int);
    
    void mutate_branch(RandomPark*, Node *&);
  
    void mutate_constant(RandomPark*, Node *&, double);
    
    void mutate_variable(RandomPark*, Node *&);
    
    void mutate_operator(RandomPark*, Node *&);
    
    void crossover(RandomPark*, Node *&, Node *&, int);
    
    void swap(Node *&, Node *&);
    
    void init(double, double, int, int);
    
    Node *copy(Node *);
    
    Node *from_buffer(char **, int&);
    
    //create tree from input expression
    Node *from_string(char*);
    void infix2postfix();
    void string2tokens();  
    int op_priority(char*);
    bool pop(char[250], char **);
    bool push(char[250], char **);
    void write_tokens(char **);
    bool stack_empty(char **);
    int stack_peek(char **);
    Node *postfix2tree();
    void clean_string(char *);
    
    State *state;
    inline void set_state(State * st){state = st;}
    
    int max_depth;
    
  private:
    double const_lo;
    
    double const_hi;
    
    double range;
    
    int num_variables;
    
    int nch, npr;
    
    char *input_line;
    char** tokens;
    char** op_stack;
    char** postfix;
    
    void clear_tokens(char**);
    
    
    
    Node **children;
    Node **parents;
    
  };
}
#endif
