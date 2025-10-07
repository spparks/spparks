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

#ifdef COMMAND_CLASS
CommandStyle(create_sites,CreateSites)

#else

#ifndef SPK_CREATE_SITES_H
#define SPK_CREATE_SITES_H

#include "pointers.h"

#ifdef SPPARKS_MAP
#include <map>
#elif SPPARKS_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

namespace SPPARKS_NS {

class CreateSites : protected Pointers {
  friend class ReadSites;

 public:
  CreateSites(class SPPARKS *);
  void command(int, char **);
  void read_sites(class AppLattice *);
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

  // geometric info for a simple, regular lattice
  // xyz me = lattice index bounds of my subdomain

  int xlo_me,xhi_me,ylo_me,yhi_me,zlo_me,zhi_me;

  // site info to pack into a message from random lattices

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

#ifdef SPPARKS_MAP
  typedef std::map<tagint,int> MyHash;
  typedef std::map<tagint,int>::iterator MyIterator;
#elif SPPARKS_UNORDERED_MAP
  typedef std::unordered_map<tagint,int> MyHash;
  typedef std::unordered_map<tagint,int>::iterator MyIterator;
#else
  typedef std::tr1::unordered_map<tagint,int> MyHash;
  typedef std::tr1::unordered_map<tagint,int>::iterator MyIterator;
#endif

  // local methods

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

/* ERROR/WARNING messages:

E: Create_sites command before app_style set

Self-explanatory.

E: Create_sites command before simulation box is defined

Self-explanatory.

E: Cannot create sites after sites already exist

Self-explanatory.

E: Cannot create sites with undefined lattice

Must use lattice commands first to define a lattice.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Create_sites region ID does not exist

Self-explanatory.

E: Creating a quantity application does not support

The application defines what variables it supports.  You cannot set a
variable with the create_sites command for a variable that isn't
supported.

E: Must use value option before basis option in create_sites command

Self-explanatory.

E: Cannot use create_sites basis with random lattice

Self-explanatory.

E: Periodic box is not a multiple of lattice spacing

UNDOCUMENTED

E: Did not create correct number of sites

One or more created sites were not assigned to a processor
correctly.

E: Bad neighbor site ID

UNDOCUMENTED

E: Random lattice has no connectivity

The cutoff distance is likely too short.

E: Ghost site was not found

Internal SPPARKS error.  Should not occur.

E: Ghost connection was not found

Internal SPPARKS error.  Should not occur.

E: Incorrect lattice neighbor count

Internal SPPARKS error.

*/
