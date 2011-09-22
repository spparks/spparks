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

#ifdef COMMAND_CLASS
CommandStyle(create_sites,CreateSites)

#else

#ifndef SPK_CREATE_SITES_H
#define SPK_CREATE_SITES_H

#include "pointers.h"

namespace SPPARKS_NS {

class CreateSites : protected Pointers {
  friend class ReadSites;

 public:
  CreateSites(class SPPARKS *);
  void command(int, char **);
  void ghosts_from_connectivity(class AppLattice *, int);

 private:
  int style,nregion,valueflag,valueindex,ivalue;
  double dvalue;
  int *basisflag,*basis_ivalue;
  double *basis_dvalue;
  int maxneigh;

  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  tagint **idneigh;            // global indices of neighbors of each site
                               // same as AppLattice neighbor, but tagint,
                               // tmp storage until convert to local indices

  int **siteijk;               // global indices of each site
                               // 0,1,2 = i,j,k lattice indices
                               // 3 = which basis atom in unit cell

  int ***cmap;                 // connectivity map for regular lattices
                               // cmap[nbasis][maxneigh][4]
                               // 0,1,2 = i,j,k lattice unit cell offsets
                               // 3 = which basis atom in unit cell


  // geometry info for building structured lattice with neighbors

  int nx,ny,nz;
  int xlo,xhi,ylo,yhi,zlo,zhi;
  double xorig,yorig,zorig;
  int latstyle,nbasis;
  double xlattice,ylattice,zlattice;

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  void structured_lattice();
  void structured_connectivity();
  void random_sites();
  void random_connectivity();

  void offsets(double **);
  void offsets_2d(int, double **, double, double, int, int **);
  void offsets_3d(int, double **, double, double, int, int **);
};

}

#endif
#endif
