/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

namespace SPPARKS_NS {

class Neighborhood{
  protected:
    int size;
    int max_nb;
    class State *state;
    int state_flag;
    char name[30];
    
  public:
    Neighborhood(char *, class State *);
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

