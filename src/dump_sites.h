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

#ifdef DUMP_CLASS

DumpStyle(sites,DumpSites)

#else

#ifndef SPK_DUMP_SITES_H
#define SPK_DUMP_SITES_H

#include "dump_text.h"

namespace SPPARKS_NS {

class DumpSites : public DumpText {
 public:
  DumpSites(class SPPARKS *, int, char **);
  ~DumpSites() {}
  void write_header(bigint, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
