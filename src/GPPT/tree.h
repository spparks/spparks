/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef TREE_H
#define TREE_H

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
    
    class Node *root;
    class Node *build_tree(class RandomPark*, int);  
    class Node *random_node(class RandomPark*, int);
  
    void mutate_branch(class RandomPark*, class Node *&);
    void mutate_constant(class RandomPark*, class Node *&, double);   
    void mutate_variable(class RandomPark*, class Node *&);   
    void mutate_operator(class RandomPark*, class Node *&);
    
    void crossover(class RandomPark*, class Node *&, class Node *&, int);
    
    void swap(class Node *&, class Node *&);
    
    void init(double, double, int, int);
    
    class Node *copy(class Node *);
    
    class Node *from_buffer(char **, int&);
    
    //create tree from input expression

    class Node *from_string(char*);
    void infix2postfix();
    void string2tokens();  
    int op_priority(char*);
    bool pop(char[250], char **);
    bool push(char[250], char **);
    void write_tokens(char **);
    bool stack_empty(char **);
    int stack_peek(char **);
    class Node *postfix2tree();
    void clean_string(char *);
    

    int state_flag;
    inline void set_state(class State * st)
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
    
    class Node **children;
    class Node **parents;

    class State *state;

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
