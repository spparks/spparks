/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
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
  double xlo, xhi, ylo, yhi, zlo, zhi ;  // bounds of cluster
  int nneigh;         // number of neighbor clusters
  int* neighlist;     // list of neighbor cluster ids
  int* pbcflags;      // list of pbc flags (-1,0,1) 
  int pbcflagsself[3];// integer pbc shifts accrued from parents

  Cluster(int, int, double, double, double, double, double, double, double, double, double, double, double, int, double*, double*);
  Cluster& operator=(const Cluster&);
  ~Cluster();
  void add_neigh(int);                     // add an entry to neighlist
  void add_pbcflags(int, int*);            // update pbc flags for neighbor
  void print(FILE*);
};

}

#endif
