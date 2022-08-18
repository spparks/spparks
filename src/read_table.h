/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
   
   Class ReadTable - added by Eric Homer, ehomer@sandia.gov
   Sep 22, 2010 - This class was largely copied from the ReadSites class.

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
CommandStyle(read_table,ReadTable)

#else

#ifndef SPK_READ_TABLE_H
#define SPK_READ_TABLE_H

#include "pointers.h"
#include "table.h"

namespace SPPARKS_NS {

class ReadTable : protected Pointers {
 public:
  ReadTable(class SPPARKS *);
  ~ReadTable();
  void command(int, char **);

 private:
  int me;
  char *line,*buffer,*keyword;
  FILE *fp;
  int compressed;
  Table *table;
  
  void header();
  void values();
  void open(char *);
  void parse_keyword(int);
  int count_words(char *);
};

}

#endif
#endif
