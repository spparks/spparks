/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef FITNESS_H
#define FITNESS_H

#include "sysptr.h"
#include "stdio.h"
#include <string>
#include "node.h"

using namespace std;
namespace SPPARKS_NS {

class Fitness : protected SysPtr {
 public:
    Fitness(SPPARKS *);
    ~Fitness();
    double compute(Node *);
    void init();
    inline double tree_value(Node *root_in){return root_in->go(var);}
    inline int get_num_var(){return nvar;}
    void set_fitness(int, char **, int);
    double get_distance();
    void update_weight(int);
    
  private:
    class RandomPark *random;

    double *var;
    int nvar;
    string fitfile_name;
    FILE *fitfile;

    void read_fit_file(int);
    int nsnaps;
    double **snaps;
    double *energy;
    int *npairs; 
    int *hits;
    double *weight;

    double *test_e;
  };
}
#endif
