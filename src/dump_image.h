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
