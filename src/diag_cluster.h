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

#ifdef DIAG_CLASS
DiagStyle(cluster,DiagCluster)

#else

#ifndef SPK_DIAG_CLUSTER_H
#define SPK_DIAG_CLUSTER_H

#include "stdio.h"
#include <stack>
#include "diag.h"

namespace SPPARKS_NS {

class DiagCluster : public Diag {
 public:
  DiagCluster(class SPPARKS *, int, char **);
  virtual ~DiagCluster();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 protected:
  int* cluster_ids;
  int ncluster, ncluster_reduced;
  double vav,rav;
  class Cluster *clustlist;
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
  int first_run;

  enum DumpStyles {STANDARD,OPENDX};
  int dump_style;
  int idump;
  char* opendxroot;
  int opendxcount;

  void analyze_clusters();
  void write_header();
  void dump_clusters(double);
  void dump_clusters_detailed();
  void generate_clusters();
  void add_cluster(int, int, double, double, int, double*);
  void free_clustlist();
};

}

#endif
#endif
