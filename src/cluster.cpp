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
		 double cx_in, double cy_in, double cz_in,
		 int nn, double* neighs, double* pbcs) {
  global_id = id;
  ivalue = iv;
  dvalue = dv;
  volume = vol;
  cx = cx_in;
  cy = cy_in;
  cz = cz_in;
  pbcflagsself[0] = pbcflagsself[1] = pbcflagsself[2] = 0;
  nneigh = nn;
  if (nneigh == 0) {
    neighlist = NULL;
    pbcflags = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    pbcflags = (int*) malloc(3*nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = static_cast<int> (neighs[i]);
      pbcflags[3*i] = static_cast<int> (pbcs[3*i]);
      pbcflags[3*i+1] = static_cast<int> (pbcs[3*i+1]);
      pbcflags[3*i+2] = static_cast<int> (pbcs[3*i+2]);
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
  cx = c.cx;
  cy = c.cy;
  cz = c.cz;
  pbcflagsself[0] = pbcflagsself[1] = pbcflagsself[2] = 0;
  nneigh = c.nneigh;
  if (nneigh == 0) {
    neighlist = NULL;
    pbcflags = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    pbcflags = (int*) malloc(3*nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = c.neighlist[i];
      pbcflags[3*i] = c.pbcflags[3*i];
      pbcflags[3*i+1] = c.pbcflags[3*i+1];
      pbcflags[3*i+2] = c.pbcflags[3*i+2];
    }
  }
  return *this;
}
  
Cluster::~Cluster() {
  if (neighlist != NULL) {
    free(neighlist);
    free(pbcflags);
  }
}

// Include cluster id in list of neighbors of this cluster
  
void Cluster::add_neigh(int id) {
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    if (id == neighlist[ineigh]) return;
  }
  nneigh++;
  neighlist = (int*) realloc(neighlist,nneigh*sizeof(int));
  pbcflags = (int*) realloc(pbcflags,3*nneigh*sizeof(int));
  neighlist[nneigh-1] = id;
  pbcflags[3*(nneigh-1)] = 0;
  pbcflags[3*(nneigh-1)+1] = 0;
  pbcflags[3*(nneigh-1)+2] = 0;
}

// Record which periodic boundaries separate 
// this cluster and the neighbor cluster.
// input id is the id of neighbor.
// pbcflags_in take values -1,0,+1.
// non-zero values are copied to pbcflags
// overwriting previous value which is
// either zero or identical to pbcflags_in.
// exception case is if pbcflags is already 
// non-zero and opposite. This means
// that the two clusters meet at both left
// and right processor boundaries, which
// can only happen for an infinite periodic cluster
// and two processors in that direction.
// centroid is not defined for infinite cluster
// so we do not need to worry about it.
//
// If it is desired to detect infinite periodic clusters,
// three different cases have to be handled,
// differing in the number of processors in
// the spanning direction, say Px.
//
// 1. Px = 1. This can only be detected
// in AppLattice::push_connected_neighbors(). If a 
// neighbor has already been assigned to the cluster
// and its distance in x is greater than lx/2, then
// the cluster has wrapped around the periodic 
// boundary in x, and so is infinite in x dir.
//
// 2. Px = 2. This can only be detected in
// Cluster::add_pbcflags(), as mentioned above.
//
// 3. Px > 3. This can only be detected at
// the level of the cluster analysis. If a 
// neighbor cluster has already been processed
// and its centroid distance is more than lx/2
// then this has "probably" wrapped around. It
// is easy to construct counter-examples. A more
// rigorous test requires performing the same test
// when the cluster is mapped onto the 1-D lattice 
// of the processor grid.

void Cluster::add_pbcflags(int id, int* pbcflags_in) {
  int ix = pbcflags_in[0];
  int iy = pbcflags_in[1];
  int iz = pbcflags_in[2];

  if (!(ix || iy || iz)) return;
  int ineigh;
  for (ineigh = 0; ineigh < nneigh; ineigh++) {
    if (id == neighlist[ineigh]) break;
  }
  if (ix) pbcflags[3*ineigh] = ix;
  if (iy) pbcflags[3*ineigh+1] = iy;
  if (iz) pbcflags[3*ineigh+2] = iz;
} 
  
void Cluster::print(FILE* fp) {
  fprintf(fp,"%d %d %g %g %g %g %g %d ",global_id,ivalue,dvalue,volume,cx,cy,cz,nneigh);
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    fprintf(fp,"%d ",neighlist[ineigh]);
  }
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    fprintf(fp,"%d %d %d ",pbcflags[3*ineigh],pbcflags[3*ineigh+1],pbcflags[3*ineigh+2]);
  }
  fprintf(fp,"\n");
}
  
