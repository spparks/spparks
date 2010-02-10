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
 public:
  CreateSites(class SPPARKS *);
  void command(int, char **);

 private:
  int style,nregion,valueflag,valueindex,ivalue;
  double dvalue;
  int *basisflag,*basis_ivalue;
  double *basis_dvalue;

  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  int latstyle,nbasis,nx,ny,nz;
  double latconstx,latconsty,latconstz;

  int ***cmap;                 // connectivity map for regular lattices
                               // cmap[nbasis][maxneigh][4]
                               // 0,1,2 = x,y,z offsets in unit cell
                               // 3 = which atom in offset unit cell

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  void structured_lattice();
  void structured_connectivity();
  void random_sites();
  void random_connectivity();
  void ghosts_from_connectivity(int);

  int connect(int, int);
  void offsets(int, double **);
  void offsets_2d(int, double **, double, double, double, double, int, int **);
  void offsets_3d(int, double **, double, double, double, double, double,
		  int, int **);
};

}

#endif
#endif
