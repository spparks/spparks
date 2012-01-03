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
CommandStyle(read_sites,ReadSites)

#else

#ifndef SPK_READ_SITES_H
#define SPK_READ_SITES_H

#include "pointers.h"

namespace SPPARKS_NS {

class ReadSites : protected Pointers {
 public:
  ReadSites(class SPPARKS *);
  ~ReadSites();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int narg,maxarg,compressed;
  char **arg;

  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appoff;

  int maxneigh;
  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;

  void open(char *);
  void header();
  void parse_keyword(int);
  void parse_coeffs(int, char *);

  void sites();
  void neighbors();
  void values();

  int count_words(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Read_sites command before app_style set

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Cannot run 2d simulation with nonperiodic Z dimension

UNDOCUMENTED

E: Cannot run 1d simulation with nonperiodic Y or Z dimension

UNDOCUMENTED

E: Cannot read Sites after sites already exist

UNDOCUMENTED

E: Cannot read Neighbors after sites already exist

UNDOCUMENTED

E: Can only read Neighbors for on-lattice applications

UNDOCUMENTED

E: Cannot read Neighbors unless max neighbors is set

UNDOCUMENTED

E: Must read Sites before Neighbors

Self-explanatory.

E: Cannot read Values before sites exist or are read

UNDOCUMENTED

E: Unknown identifier in data file: %s

Self-explanatory.

E: Site file has no Sites, Neighbors, or Values

UNDOCUMENTED

E: No Sites defined in site file

UNDOCUMENTED

E: No Neighbors defined in site file

UNDOCUMENTED

E: Unexpected end of data file

Self-explanatory.

E: Data file dimension does not match existing box

UNDOCUMENTED

E: Data file number of sites does not match existing sites

UNDOCUMENTED

E: Off-lattice application data file cannot have maxneigh setting

UNDOCUMENTED

E: Data file maxneigh setting does not match existing sites

UNDOCUMENTED

E: Data file simluation box different that current box

UNDOCUMENTED

E: System in site file is too big

UNDOCUMENTED

E: Incorrect site format in data file

Self-explanatory.

E: Did not assign all sites correctly

One or more sites in the read_sites file were not assigned to 
a processor correctly.

E: Invalid site ID in Sites section of data file

Self-explanatory.

E: Too many neighbors per site

Internal SPPARKS error.

E: Incorrect value format in data file

Self-explanatory.

E: Cannot open gzipped file

Self-explantory.

E: Cannot open file %s

Self-explanatory.

*/
