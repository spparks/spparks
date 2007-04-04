/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef TREE_H
#define TREE_H

#include "random_park.h"
#include "node.h"
#include <cmath>
#include "sysptr.h"

namespace SPPARKS {

#define MUT_CONST    0
#define MUT_VAR      1
#define MUT_OPER     2
#define MUT_BRANCH   3
#define MUT_CROSS    4
#define MUT_SWAP     5

class Tree : protected SysPtr {
 public:

    Tree();
    ~Tree();
    
    Node *root;
    
    Node *build_tree(int);
    
    Node *random_node(int);

    void mutate_branch(Node *&);

    void mutate_constant(Node *&, double);

    void mutate_variable(Node *&);

    void mutate_operator(Node *&);

    void crossover(Node *&, Node *&, int);

    void swap(Node *&, Node *&);

    void init(int, double, double, int, int);

    Node *copy(Node *);

    Node *from_buffer(char **, int&);

    int max_depth;

  private:
    class RandomPark *random;
    
    double const_lo;
    
    double const_hi;
    
    double range;
    
    int num_variables;

    int nch, npr;

    Node **children;
    Node **parents;

  };
}
#endif
