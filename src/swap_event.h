/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SWAP_EVENT_H
#define SWAP_EVENT_H

#include "event.h"

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

