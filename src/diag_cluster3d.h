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

#ifndef DIAG_CLUSTER3D_H
#define DIAG_CLUSTER3D_H

#include "stdio.h"
#include <stack>
#include "cluster.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagCluster3d : public Diag {
  friend class SweepLattice3d;

 public:
  DiagCluster3d(class SPPARKS *, int, char **);
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

  int nx_global,ny_global,nz_global,nx_procs,ny_procs,nz_procs,delpropensity;
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
