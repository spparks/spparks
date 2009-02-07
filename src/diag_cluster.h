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

#ifndef DIAG_CLUSTER_H
#define DIAG_CLUSTER_H

#include "stdio.h"
#include <stack>
#include "cluster.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagCluster : public Diag {
 public:
  DiagCluster(class SPPARKS *, int, char **);
  virtual ~DiagCluster();

  void init(double);
  void compute(double, int, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  // Functions and Data for Cluster Analysis
  void analyze_clusters(double);
  void write_header();
  void dump_clusters(double);
  void dump_clusters_detailed();
  void generate_clusters();
  void add_cluster(int, int, double, double, int, double*);
  void free_clustlist();

  int* cluster_ids;
  int ncluster, ncluster_reduced;
  Cluster* clustlist;
  std::stack<int> cluststack;      // stack for performing cluster analysis

  FILE *fp, *fpdump;
  class AppLattice *applattice;

  class CommLattice *comm;

  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // simulation box bounds

  int *id;                     // global ID (1-N) of site
  double **xyz;                // coords of site

  int nglobal;                 // global # of sites
  int nlocal;                  // # of sites I own
  int nghost;                  // # of ghost sites I store

  int nx_global,ny_global,nz_global;

  enum DumpStyles {STANDARD,OPENDX};
  int dump_style;
  int idump;
  char* opendxroot;
  int opendxcount;
};

}

#endif
