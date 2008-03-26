
#ifndef EVENT_H
#define EVENT_H

#include "random_park.h"
#include "state.h"
#include <string>
#include "sysptr.h"
//#include "expression.h"
//#include "neighborhood.h"
#include <iostream>

using namespace std;
namespace SPPARKS {

#define FLIP        0
#define KFLIP       1
#define SWAP        2
#define KRING       3

#define SRC_TO_NBR  10
#define SRC_ENUM    11

#define PROP_BOLTZ  0
#define PROP_EVAL   1
#define PROP_NB_SUM 2
  
  class Event{
    
  public:
    Event();
    virtual ~Event();
    int type;
    
    inline char* get_name(){return &name[0];}
    
    inline void set_name(char *nm_in)
      {strcpy(name,nm_in);}
    inline void set_variable(int type, int *ip)
      {vtype = type; ivar=ip;}
    inline void set_variable(int type, double *dp)
      {vtype = type; dvar=dp;}
    inline void set_neighborhood(class Neighborhood * nbr_in)
      {nbr = nbr_in;}

    inline void set_enum_bounds(int ilo, int ihi)
      {lo = ilo; hi = ihi;}
    inline void set_src(int srctype_in){srctype = srctype_in;}
    virtual int site_event(int, double)=0;
    void set_propensity(int, class Xpression *);
    inline void set_temperature(double t)
      {if(t>0)beta = 1.0/t;else beta = 0.0;}
    virtual void init(){}

    virtual double site_propensity(int site){return 0;}

    class Neighborhood *nbr;

    
  protected:
    
    char name[30];

    int vtype;
    double *dvar;
    int *ivar;
    int hi, lo;
    int srctype;
    int maxnb;

    class Xpression *prop;
    int ptype;
    double beta;
    
  };
}

#endif

