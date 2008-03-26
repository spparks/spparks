#include "swap_event.h"
#include "event.h"
#include <iostream>
#include "random_park.h"
#include "state.h"
#include "neighborhood.h"
#include "expression.h"

using namespace std;
using namespace SPPARKS;


/* ---------------------------------------------------------------------- */
SwapEvent::SwapEvent() : Event()
{
  type = SWAP;
  ptype = PROP_EVAL;
}
/* ---------------------------------------------------------------------- */
SwapEvent::~SwapEvent() 
{
  delete [] endval;
}

/* ---------------------------------------------------------------------- */
void SwapEvent::init()
{
  
}
/* ---------------------------------------------------------------------- */
int SwapEvent::site_event(int site, double threshold)
{
  double tot_prop = 0;
  int val = -1000;
  int w;
  double eold = 0;
  double enew = 0;


  if (srctype == SRC_TO_NBR){
    int nnb = nbr->numneigh[site];
    for (int v = 0; v < nnb; v++){
      int nb = nbr->neighbor[site][v];
      if (ptype == PROP_BOLTZ) 
	eold = prop->sum(site)+prop->sum(nb);

      val = ivar[nb];
      ivar[nb] = ivar[site];
      ivar[site] = val;

      if     (ptype == PROP_BOLTZ)  {
	double enew = prop->sum(site) + prop->sum(nb);
	if (enew <= eold) tot_prop += 1.0;
	else tot_prop += exp((eold-enew)*beta);
      }
      else if(ptype == PROP_NB_SUM) 
	tot_prop += prop->sum(site) + prop->sum(nb);
      else if(ptype == PROP_EVAL)   
	tot_prop += prop->eval(site) + prop->eval(nb);
      
      if (tot_prop>= threshold) return 100 + nb; //place of swap nbr (-100)
      
      val = ivar[nb];
      ivar[nb] = ivar[site];
      ivar[site] = val;
    }
  }
  
  
  return 0;
}
/* ---------------------------------------------------------------------- */
double SwapEvent::site_propensity(int site)
{
  double tot_prop = 0;
  int val = -1;
  int w;
  double eold;

  if (srctype == SRC_TO_NBR){
    int nnb = nbr->numneigh[site];
    for (int v = 0; v < nnb; v++){
      int nb = nbr->neighbor[site][v];
      if (ptype == PROP_BOLTZ) 
	eold = prop->sum(site)+prop->sum(nb);

      val = ivar[nb];
      ivar[nb] = ivar[site];
      ivar[site] = val;

      if     (ptype == PROP_BOLTZ)  {
	double enew = prop->sum(site) + prop->sum(nb);
	if (enew <= eold) tot_prop += 1.0;
	else tot_prop += exp((eold-enew)*beta);
      }
      else if(ptype == PROP_NB_SUM) tot_prop += prop->sum(site);
      else if(ptype == PROP_EVAL)   tot_prop += prop->eval(site);

      val = ivar[nb];
      ivar[nb] = ivar[site];
      ivar[site] = val;
    }
  }

  return tot_prop;
}
/* ---------------------------------------------------------------------- */
