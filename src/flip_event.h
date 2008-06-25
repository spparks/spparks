/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef FLIP_EVENT_H
#define FLIP_EVENT_H

#include "event.h"

namespace SPPARKS_NS {
  
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

