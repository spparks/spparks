
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
    std::string nbname;
    
  public:
    Neighborhood(char *, State *);
    ~Neighborhood();
    
    int **nbr;
    
    void init(char *);
  };
}

#endif

