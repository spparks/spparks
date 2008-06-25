/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

/*-----------------------------------------------------------------------
An event requires a propensity expression and is an update
method. The propensity expression is supplied. 
The event is referenced to a particular lattice site in this context.
A sweep of sites visits each site in turn and calculates the propensities
of the individual events.
An event is defined by the local update of the state that it performs:
1. Flip -- change the value of an attribute on-site.
2. K-flip -- change the value of an attribute on k local sites
3. Swap -- interchage the values of attribute on two sites.
4. K-Ring -- permute the values of k sites  
-----------------------------------------------------------------------*/

#include "event.h"
#include "expression.h"
#include "neighborhood.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

Event::Event()
{
  hi = 0;
  lo = 0;
  type = -1;
  vtype = 0;
  ptype = 0;
  nbr = NULL;
  dvar = NULL;
  ivar = NULL;
}

/* ---------------------------------------------------------------------- */

Event::~Event()
{
  delete nbr;
}

/* ---------------------------------------------------------------------- */

void Event::set_propensity(int ptp, class Xpression *xpr)
{
  ptype = ptp; 
  prop = xpr;  
  if(nbr!=NULL) prop->set_neighborhood(nbr);
}
