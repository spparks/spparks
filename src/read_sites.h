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
