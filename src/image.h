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

#ifndef SPK_IMAGE_H
#define SPK_IMAGE_H

#include "math.h"
#include "stdio.h"
#include "pointers.h"

namespace SPPARKS_NS {

class Image : protected Pointers {
 public:
  int width,height;             // size of image
  double theta,phi;             // view image from theta,phi
  double xctr,yctr,zctr;        // center of image in user coords
  double up[3];                 // up direction in image
  double zoom;                  // zoom factor
  double persp;                 // perspective factor
  double shiny;                 // shininess of objects
  int ssao;                     // SSAO on or off
  int seed;                     // RN seed for SSAO
  double ssaoint;               // strength of shading from 0 to 1
  double *boxcolor;             // color to draw box outline with
  int background[3];            // RGB values of background

  Image(class SPPARKS *);
  ~Image();
  void buffers();
  void clear();
  void merge();
  void write_JPG(FILE *);
  void write_PPM(FILE *);
  void view_params(double, double, double, double, double, double);

  void color_minmax(int, double *, int);
  void draw_sphere(double *, double *, double);
  void draw_cube(double *, double *, double);
  void draw_cylinder(double *, double *, double *, double, int);
  void draw_triangle(double *, double *, double *, double *);
  void draw_box(double (*)[3], double);
  void draw_axes(double (*)[3], double);

  int colormap(int, char **);
  int addcolor(char *, double, double, double);
  double *element2color(char *);
  double element2diam(char *);
  double *value2color(double);
  double *color2rgb(const char *, int index=0);
  int default_colors();

 private:
  int me,nprocs;
  int npixels;

  double *depthBuffer,*surfaceBuffer;
  double *depthcopy,*surfacecopy;
  char *imageBuffer,*rgbcopy,*writeBuffer;

  // constant view params

  double FOV;
  double ambientColor[3];

  double keyLightTheta;
  double keyLightPhi;
  double keyLightColor[3];

  double fillLightTheta;
  double fillLightPhi;
  double fillLightColor[3];

  double backLightTheta;
  double backLightPhi;
  double backLightColor[3];

  double specularHardness;
  double specularIntensity;

  double SSAORadius;
  int SSAOSamples;
  double SSAOJitter;

  // dynamic view params

  double zdist;
  double tanPerPixel;
  double camDir[3],camUp[3],camRight[4],camPos[3];
  double keyLightDir[3],fillLightDir[3],backLightDir[3];
  double keyHalfDir[3];

  // color values

  int ncolors;
  char **username;
  double **userrgb;

  // color map

  int mstyle,mrange;               // 2-letter style/range of color map
  int mlo,mhi;                     // bounds = NUMERIC or MINVALUE or MAXVALUE
  double mlovalue,mhivalue;        // user bounds if NUMERIC
  double locurrent,hicurrent;      // current bounds for this snapshot
  double mbinsize,mbinsizeinv;     // bin size for sequential color map

  struct MapEntry {
    int single,lo,hi;              // NUMERIC or MINVALUE or MAXVALUE
    double svalue,lvalue,hvalue;   // actual value
    double *color;                 // RGB values
  };

  MapEntry *mentry;
  int nentry;
  double interpolate[3];

  // SSAO RNG

  class RandomPark *random;

  // internal methods

  void draw_pixel(int, int, double, double *, double*);
  void compute_SSAO();

  // inline functions

  inline double saturate(double v) {
    if (v < 0.0) return 0.0;
    else if (v > 1.0) return 1.0;
    else return v;
  }

  inline double distance(double* a, double* b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) + 
		(a[1] - b[1]) * (a[1] - b[1]) + 
		(a[2] - b[2]) * (a[2] - b[2]));
  }
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid image up vector

UNDOCUMENTED

E: Invalid image color range

UNDOCUMENTED

*/
