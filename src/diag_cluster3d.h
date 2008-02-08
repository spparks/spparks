/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_CLUSTER3D_H
#define DIAG_CLUSTER3D_H

#include "stdio.h"
#include <stack>
#include "cluster.h"
#include "diag.h"

namespace SPPARKS {

class DiagCluster3d : public Diag {
  friend class SweepLattice3d;

 public:
  DiagCluster3d(class SPK *, int, char **);
  virtual ~DiagCluster3d();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  // Functions and Data for Cluster Analysis
  void analyze_clusters(double);
  void write_header();
  void dump_clusters(double);
  void dump_clusters_detailed();
  void generate_clusters();
  void add_cluster(int, double, int, double*);
  void free_clustlist();

  int*** cluster_ids;
  int ncluster, ncluster_reduced;
  double radius;
  Cluster* clustlist;
  std::stack<int> cluststack;      // stack for performing cluster analysis

  FILE *fp, *fpdump;
  class AppLattice3d *applattice3d;

  int nx_global,ny_global,nz_global,nx_procs,ny_procs,nz_procs,delghost;
  int nx_local,ny_local,nz_local,nx_offset,ny_offset,nz_offset;
  int nxlo,nylo,nzlo,nxhi,nyhi,nzhi;

  enum DumpStyles {STANDARD,OPENDX};
  int dump_style;
  int idump;
  char* opendxroot;
  int opendxcount;
};

}

#endif
