/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "flip_event.h"
#include "neighborhood.h"
#include "expression.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

FlipEvent::FlipEvent() : Event()
{
  type = FLIP;
  ptype = PROP_EVAL;
}

/* ---------------------------------------------------------------------- */

FlipEvent::~FlipEvent() 
{
  delete [] endval;
}

/* ---------------------------------------------------------------------- */

void FlipEvent::init()
{
  if (srctype == SRC_TO_NBR){
    maxnb = nbr->maxnb;
    endval = new int[maxnb];
  }
  else if (srctype == SRC_ENUM){
    nevent = hi - lo + 1;
    endval = new int[nevent];
  }

  if (srctype == SRC_ENUM){
    for (int v = lo; v < hi; v++) endval[v-lo] = v;
  }
  
}

/* ---------------------------------------------------------------------- */

int FlipEvent::site_event(int site, double threshold)
{
  double tot_prop = 0;
  nevent = 0;
  int val = -1;
  int w;
  double eold = 0;
  double enew = 0;

  int oldval = ivar[site];
  if (ptype == PROP_BOLTZ) eold = prop->sum(site);

  if (srctype == SRC_TO_NBR) {
    int nnb = nbr->numneigh[site];

    // cout << "Neighbor values: ";
    // for (int n = 0; n < nbr->numneigh[site]; n++) 
      // cout<<" "<<ivar[nbr->neighbor[site][n]];
    // cout <<endl;
    

    for (int v = 0; v < nnb; v++){
      val = ivar[nbr->neighbor[site][v]];
      if (val == ivar[site]) continue;
      for(w = 0; w < nevent; w++) if(val == endval[w]) break;
      if (w < nevent) continue;
      endval[nevent] = val;
      nevent++; 
      ivar[site] = endval[w];
      if (ptype == PROP_BOLTZ)  {
	enew = prop->sum(site);
	if (enew <= eold) tot_prop += 1.0;
	else tot_prop += exp((eold-enew)*beta);
      }
      else if(ptype == PROP_NB_SUM) tot_prop += prop->sum(site);
      else if(ptype == PROP_EVAL)   tot_prop += prop->eval(site);

      if (tot_prop>= threshold) return 1;
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

double FlipEvent::site_propensity(int site)
{
  double tot_prop = 0;
  nevent = 0;

  int val = -1;
  int w;
  double eold;

  if (srctype == SRC_TO_NBR){
    int nnb = nbr->numneigh[site];
    for (int v = 0; v < nnb; v++){
      val = ivar[nbr->neighbor[site][v]];
      if (val == ivar[site]) continue;
      for(w=0;w<nevent;w++) if(val==endval[w])break;
      if(w<nevent)continue;
      endval[nevent] = val;
      nevent++;
    }
  }

  int oldval = ivar[site];
  if (ptype == PROP_BOLTZ) eold = prop->sum(site);

  for (int v = 0; v < nevent; v++){

    ivar[site] = endval[v];

    if (ptype == PROP_BOLTZ) {
      double enew = prop->sum(site);
      if (enew <= eold) tot_prop += 1.0;
      else tot_prop += exp((eold-enew)*beta);
    }
    else if(ptype == PROP_NB_SUM) tot_prop += prop->sum(site);
    else if(ptype == PROP_EVAL)   tot_prop += prop->eval(site);
  }

  ivar[site] = oldval;
  return tot_prop;
}
