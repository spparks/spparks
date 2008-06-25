/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef CLUSTER_H
#define CLUSTER_H

namespace SPPARKS_NS {

class Cluster {
 private:
  double volume;
  int global_id;
  int nneigh;
  int* neighlist;
  
 public:
  Cluster(int, double, int, double*);
  Cluster& operator=(const Cluster&);
  ~Cluster();
  void add_neigh(int);
  void print(FILE*);

  friend class DiagCluster;
  friend class DiagCluster2d;
  friend class DiagCluster3d;
};

}

#endif
