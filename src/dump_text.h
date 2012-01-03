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

DumpStyle(text,DumpText)

#else

#ifndef SPK_DUMP_TEXT_H
#define SPK_DUMP_TEXT_H

#include "dump.h"

namespace SPPARKS_NS {

class DumpText : public Dump {
 public:
  DumpText(class SPPARKS *, int, char **);
  virtual ~DumpText();

 protected:
  int ioptional;             // where optional trailing args start

  int *vtype;                // type of each vector (INT, DOUBLE)
  int *vindex;               // index into int,double packs
  char **vformat;            // format string for each vector element

  int iregion;               // -1 if no region, else which region
  char *idregion;            // region ID

  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshold on for each nthresh
  int *thresh_op;            // threshold operation for each nthresh
  double *thresh_value;      // threshold value for each nthresh
  int *thresh_index;         // N index for iN and dN thresholds

  char *columns;             // text describing columns of dump output

  int nchoose;               // # of selected atoms
  int maxlocal;              // size of atom selection and variable arrays
  int *choose;               // lists of sites chosen for output
  double *dchoose;           // value for each atom to threshhold against
  int *clist;                // compressed list of indices of selected atoms

  // private methods

  virtual void init_style();
  int count();
  void pack();
  void write_header(int, double);
  void write_data(int, double *);
  int parse_fields(int, char **);
  virtual int modify_param(int, char **);

  typedef void (DumpText::*FnPtrHeader)(int, double);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int, double);
  void header_text(int, double);

  typedef void (DumpText::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_text(int, double *);

  typedef void (DumpText::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_site(int);
  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_energy(int);
  void pack_propensity(int);
  void pack_iarray(int);
  void pack_darray(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid attribute in dump text command

UNDOCUMENTED

E: Dump requires propensity but no KMC solve performed

Only KMC solvers compute propensity for sites.

E: Region ID for dump text does not exist

UNDOCUMENTED

E: Dumping a quantity application does not support

The application defines what variables it supports.  You cannot
output a variable in a dump that isn't supported.

E: Invalid keyword in dump command

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Dump_modify region ID does not exist

UNDOCUMENTED

E: Threshold for a quantity application does not support

The application defines what variables it supports.  You cannot do a
threshold test with the dump command on a variable that isn't
supported.

E: Invalid dump_modify threshold operator

Self-explanatory.

*/
