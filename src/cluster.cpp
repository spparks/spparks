/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "cluster.h"

using namespace SPPARKS_NS;

Cluster::Cluster(int id, int iv, double dv, double vol,
		 int nn, double* neighs) {
  global_id = id;
  ivalue = iv;
  dvalue = dv;
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
  ivalue = c.ivalue;
  dvalue = c.dvalue;
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
  fprintf(fp,"%d %d %g %g %d ",global_id,ivalue,dvalue,volume,nneigh);
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    fprintf(fp,"%d ",neighlist[ineigh]);
  }
  fprintf(fp,"\n");
}
  
