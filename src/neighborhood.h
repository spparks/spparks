
#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "random_park.h"
#include "state.h"
#include <string>
#include "sysptr.h"

namespace SPPARKS {
  class Neighborhood{
    
  protected:
    
    int size;
    int max_nb;
    State *state;
    int state_flag;
    char name[30];
    
  public:
    Neighborhood(char *, State *);
    Neighborhood(char *, int *, int **, int, int);
    ~Neighborhood();
    
    int *numneigh_in;
    int **neighbor_in;
    int **neighbor;
    int *numneigh;
    int nsites;
    int maxnb;

    void init(int, char**);
    inline char* get_name(){return &name[0];}
  };
}

#endif

