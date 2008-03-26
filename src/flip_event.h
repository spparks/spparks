/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef FLIP_EVENT_H
#define FLIP_EVENT_H

#include "random_park.h"
#include "sysptr.h"
#include "event.h"
//#include "expression.h"
#include <iostream>

using namespace std;
namespace SPPARKS {
  
  class FlipEvent : public Event {
  public: 
    FlipEvent();
    ~FlipEvent();
    
    void init();
    int site_event(int, double);
    double site_propensity(int);

  private:
    int *endval;
    int nevent;
    
  };
}

#endif

