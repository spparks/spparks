/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef CLUSTER_H
#define CLUSTER_H

namespace SPPARKS {

class Cluster {
 private:
  int volume;
  int global_id;
  int nneigh;
  int* neighlist;
  
 public:
  Cluster(int, int, int, int*);
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
