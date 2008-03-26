/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#ifndef SWAP_EVENT_H
#define SWAP_EVENT_H

#include "random_park.h"
#include "sysptr.h"
#include "event.h"
//#include "expression.h"
#include <iostream>

using namespace std;
namespace SPPARKS {
  
  class SwapEvent : public Event {
  public: 
    SwapEvent();
    ~SwapEvent();
    
    void init();
    int site_event(int, double);
    double site_propensity(int);

  private:
    int *endval;
    int nevent;
    
  };
}

#endif

