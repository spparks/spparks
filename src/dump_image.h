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

DumpStyle(image,DumpImage)

#else

#ifndef SPK_DUMP_IMAGE_H
#define SPK_DUMP_IMAGE_H

#include "math.h"
#include "dump_text.h"

namespace SPPARKS_NS {

class DumpImage : public DumpText {
 public:
  DumpImage(class SPPARKS *, int, char**);
  ~DumpImage();

 private:
  int filetype;
  int shape,boundflag,scolor,sdiam,boundvalue,boundindex;
  int crange,drange;
  int clo,chi,dlo,dhi;
  double sdiamvalue;
  double bounddiam;

  char *thetastr,*phistr;          // variables for view theta,phi
  int thetavar,phivar;             // index to theta,phi vars
  int cflag;                       // static/dynamic box center
  double cx,cy,cz;                 // fractional box center
  char *cxstr,*cystr,*czstr;       // variables for box center
  int cxvar,cyvar,czvar;           // index to box center vars
  char *upxstr,*upystr,*upzstr;    // view up vector variables
  int upxvar,upyvar,upzvar;        // index to up vector vars
  char *zoomstr,*perspstr;         // view zoom and perspective variables
  int zoomvar,perspvar;            // index to zoom,persp vars
  int boxflag,axesflag;            // 0/1 for draw box and axes
  double boxdiam,axeslen,axesdiam; // params for drawing box and axes

  int viewflag;                    // overall view is static or dynamic

  class AppLattice *applattice;

  double *diamattribute;
  double **colorattribute;
  int *color_memflag;
  double *boundcolor;

  class Image *image;              // class that renders each image

  void init_style();
  int modify_param(int, char **);
  void write(double);

  void box_center();
  void view_params();
  void box_bounds();

  void create_image();

  void bounds(char *, int, int, int &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump image filename

UNDOCUMENTED

E: Cannot dump JPG file

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Dump image with quantity application does not support

UNDOCUMENTED

E: Invalid dump image theta value

UNDOCUMENTED

E: Dump image persp option is not yet supported

UNDOCUMENTED

W: Using dump image boundary with spheres

UNDOCUMENTED

E: Dump image boundary requires lattice app

UNDOCUMENTED

E: Dump image drange must be set

UNDOCUMENTED

E: Dump image crange must be set

UNDOCUMENTED

E: Dump image requires one snapshot per file

UNDOCUMENTED

E: Variable name for dump image theta does not exist

UNDOCUMENTED

E: Variable for dump image theta is invalid style

UNDOCUMENTED

E: Variable name for dump image phi does not exist

UNDOCUMENTED

E: Variable for dump image phi is invalid style

UNDOCUMENTED

E: Variable name for dump image center does not exist

UNDOCUMENTED

E: Variable for dump image center is invalid style

UNDOCUMENTED

E: Variable name for dump image zoom does not exist

UNDOCUMENTED

E: Variable for dump image zoom is invalid style

UNDOCUMENTED

E: Variable name for dump image persp does not exist

UNDOCUMENTED

E: Variable for dump image persp is invalid style

UNDOCUMENTED

E: Invalid dump image zoom value

UNDOCUMENTED

E: Invalid dump image persp value

UNDOCUMENTED

E: Invalid color in dump_modify command

UNDOCUMENTED

E: Dump_modify scolor requires integer attribute for dump image color

UNDOCUMENTED

E: Dump_modify sdiam requires integer attribute for dump image diameter

UNDOCUMENTED

*/
