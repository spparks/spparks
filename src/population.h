/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef POPULATION_H
#define POPULATION_H

#include "random_park.h"
#include "node.h"
#include "tree.h"
#include "sysptr.h"
#include "fitness.h"

namespace SPPARKS {

class Population : protected SysPtr {
 public:

    Population();
    ~Population();
    void init(int, int, Tree *, Fitness *, int, double);

    void build_new_population();
    int linear_select_mutation();
    bool accept(Node *, Node *, double &);

    inline double get_temperature(){return temperature;}
    inline void set_temperature(double t){
      temperature = t; 
      beta = 1./temperature;
    }

    Node*get_tree(int);
    Node*get_best_tree();
    void write_best_tree(FILE *);
    void calc_best_fit();
    inline double get_best_fitness(){return best_fit;}
    double get_stats(double &); //returns mean (stdv as argument)

    void update_weights();
    inline double *get_weights(int &n){n=MUT_SWAP+1;return weight;}

    void swap_nbr(int, double, Node*);
  private:

    void set_weights();

    RandomPark *random;
    Fitness *fit;
    Tree *tree;

    Node **roots;
    int ntrees;

    int best_index;
    double best_fit;
    double *fitness;

    double weight[MUT_SWAP+1];
    int accepted[MUT_SWAP+2];

    double temperature, beta;

    int swap_up;
    int swap_down;
  };
}
#endif
