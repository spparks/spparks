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

#ifdef DUMP_CLASS

DumpStyle(stitch, DumpStitch)

#else

#ifndef SPK_DUMP_STITCH_H
#define SPK_DUMP_STITCH_H

#include "dump.h"
#include "stitch.h"

namespace SPPARKS_NS {

class DumpStitch : public Dump {
public:
  DumpStitch(class SPPARKS *, int, char **);
  virtual ~DumpStitch();
  void write(double);

private:
  int *fields;               // descriptor for each field
  int *vtype;                // type of each vector (INT, DOUBLE)
  int *vindex;               // index into int,double columns
  int64_t *stitch_field_ids; // stitch field ids

  class AppLattice *applattice;
  StitchFile *fid;

  // private methods

  void init_style();
  int count() { return 0; } // NOT USED by 'dump_stitch'
  void pack(tagint *) {}
  virtual void write_header(bigint, double) {}
  virtual void write_data(int, double *) {}
  int parse_fields(int, char **);
  void create_stitch_field_ids();
};

} // namespace SPPARKS_NS

#endif
#endif

    /* ERROR/WARNING messages:

    */
