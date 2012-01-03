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
  int *cluster_ids;
  int ncluster,ncluster_reduced;
  double vav,rav;
  class Cluster *clustlist;
  std::stack<int> cluststack;      // stack for performing cluster analysis

  FILE *fp, *fpdump;
  class AppLattice *applattice;

  class CommLattice *comm;

  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // simulation box bounds

  double **xyz;                // coords of site
  tagint *idsite;              // global site ids

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

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

E: Cannot open diag_style cluster output file

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot open diag_style cluster dump file

Self-explanatory.

E: Cannot use diag_style cluster without a lattice defined

This diagnostic uses the lattice style to dump OpenDx files.

E: Diag_style cluster incompatible with lattice style

UNDOCUMENTED

E: Diag_style cluster nx,ny,nz = 0

UNDOCUMENTED

E: Diag cluster does not work if ncluster > 2^31

UNDOCUMENTED

E: Mismatch in counting for dbufclust

Self-explanatory.

E: Diag cluster ivalue in neighboring clusters do not match

Internal SPPARKS error.

E: Diag cluster dvalue in neighboring clusters do not match

Internal SPPARKS error.

E: Diag style cluster dump file name too long

Self-explanatory.

E: Cannot open diag style cluster dump file

Self-explanatory.

*/
