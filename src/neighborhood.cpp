#include "neighborhood.h"
#include <iostream>
#include "random_park.h"
#include "state.h"


using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */
Neighborhood::Neighborhood(char *name_in, State *state_in)
{
  memset(name,0,30);
  strcpy(name,name_in);
  state = state_in;
  state_flag = 1;
  size = state->get_size();
  maxnb = 0;
}
/* ---------------------------------------------------------------------- */
Neighborhood::Neighborhood(char *name_in, int *num, 
			   int **nghb, int ns, int maxn_in)
{
  memset(name,0,30);
  strcpy(name,name_in);
  state_flag = 0;

  numneigh_in = num;
  neighbor_in = nghb;
  nsites = ns;
  maxnb = maxn_in;
}
/* ---------------------------------------------------------------------- */
Neighborhood::~Neighborhood()
{

}
/* ---------------------------------------------------------------------- */
void Neighborhood::init(int narg, char **arg)
{
  int order = 0;
  
  //parse description of neighborhood
  if(state_flag == 0){      // lattice_app
    
    if(strcmp(arg[1],"nearest")==0){
      order = atoi(arg[2]);
      if(order == 1){   //copy the neighbor structure from the lattice
	neighbor = new int*[nsites];
	numneigh = new int[nsites];
	for(int s = 0; s < nsites; s++){
	  numneigh[s] = numneigh_in[s];
	  neighbor[s] = new int[numneigh[s]];
	  for(int n = 0; n < numneigh[s]; n++)
	    neighbor[s][n] = neighbor_in[s][n];
	  
	}
      }
    }
  }
  
  else if (state_flag == 1){     // template app
    
  }
  
}
  /* ---------------------------------------------------------------------- */
