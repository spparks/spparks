/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "cluster.h"

using namespace SPPARKS;

Cluster::Cluster(int id, double vol, int nn, double* neighs) {
  global_id = id;
  volume = vol;
  nneigh = nn;
  if (nneigh == 0) {
    neighlist = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = static_cast<int> (neighs[i]);
    }
  }
}

// Define assigment operator.
// Because memory is allocated before
// the object is initialized, can not use constructor.
Cluster& Cluster::operator=(const Cluster& c) {
  global_id = c.global_id;
  volume = c.volume;
  nneigh = c.nneigh;
  if (nneigh == 0) {
    neighlist = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = c.neighlist[i];
    }
  }
  return *this;
}
  
Cluster::~Cluster() {
  if (neighlist != NULL) free(neighlist);
}
  
void Cluster::add_neigh(int id) {
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    if (id == neighlist[ineigh]) return;
  }
  nneigh++;
  neighlist = (int*) realloc(neighlist,nneigh*sizeof(int));
  neighlist[nneigh-1] = id;
}
  
void Cluster::print(FILE* fp) {
  fprintf(fp,"%d %g %d ",global_id,volume,nneigh);
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    fprintf(fp,"%d ",neighlist[ineigh]);
  }
  fprintf(fp,"\n");

}
  
