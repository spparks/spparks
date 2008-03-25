/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef TREE_H
#define TREE_H

#include "node.h"
#include <cmath>
#include "random_park.h"
#include "state.h"
#include "sysptr.h"
#include <iostream>

using namespace std;
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
    

    int state_flag;
    inline void set_state(State * st)
      {state_flag == 1; state = st;}

    int lattice_flag;
    inline void set_lattice(int **iarr, double **darr)
      {
	if(lattice_flag == 0){
	  lattice_flag = 1; 
	  iarray = iarr;
	  darray = darr;
	}
      }   
    int get_variable(char *, int &);
    inline void set_names(char **intnm, int ni, char **dblnm, int di)
      {
	name_int = intnm; n_int_var = ni;
	name_dbl = dblnm; n_dbl_var = di;
      }

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

    State *state;
    //lattice vars 
    double **darray;
    int **iarray;
    char **name_int;
    char **name_dbl;
    int n_int_var;
    int n_dbl_var;

  };
}
#endif
