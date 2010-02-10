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
  int global_id;
  int ivalue;
  double dvalue;
  double volume;
  int nneigh;
  int* neighlist;

  Cluster(int, int, double, double, int, double*);
  Cluster& operator=(const Cluster&);
  ~Cluster();
  void add_neigh(int);
  void print(FILE*);
};

}

#endif
