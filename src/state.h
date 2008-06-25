/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef STATE_H
#define STATE_H

#define NVARS_MAX 100

#define BOOL_TYPE   0
#define INT_TYPE    1
#define DOUBLE_TYPE 2

#include "random_park.h"
#include "sysptr.h"

namespace SPPARKS_NS {
  
class State {
  protected:
    int nvars_current;
    int size;
    int step[NVARS_MAX];
    int data_type[NVARS_MAX];
    char name[NVARS_MAX][16];
    void *attr[NVARS_MAX];
    
    RandomPark *rnd;
    int seed;
    
  public:
    State(int, int);
    ~State();
    
    void add_attribute(char *, int);
    void *get_attribute(char *, int&);
    int get_attribute(char *);
    inline char *get_name(int n){return name[n];}
    inline int get_type(int n){return data_type[n];}
    inline int get_size(){return size;}
    void parse_variable(char *);
    void init_var(int, char **);
  };
  
}

#endif

