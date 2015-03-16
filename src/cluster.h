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

#ifndef SPK_CLUSTER_H
#define SPK_CLUSTER_H

namespace SPPARKS_NS {

class Cluster {
 public:
  int global_id;      // cluster id
  int ivalue;         // integer value of cluster e.g. spin 
  double dvalue;      // floating point value of cluster e.g. 0.0
  double volume;      // number of sites in cluster
  double cx, cy, cz;  // centroid of cluster
  int nneigh;         // number of neighbor clusters
  int* neighlist;     // list of neighbor cluster ids
  int* pbcflags;      // list of pbc flags (-1,0,1) 
  int pbcflagsself[3];// integer pbc shifts accrued from parents

  Cluster(int, int, double, double, double, double, double, int, double*, double*);
  Cluster& operator=(const Cluster&);
  ~Cluster();
  void add_neigh(int);                     // add an entry to neighlist
  void add_pbcflags(int, int*);            // update pbc flags for neighbor
  void print(FILE*);
};

}

#endif
