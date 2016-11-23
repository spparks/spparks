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

#include "string.h"
#include "stdlib.h"
#include "dump_vtk.h"
#include "app.h"
#include "domain.h"
#include "lattice.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // in other dump files

enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
     FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};  // in other files

/* ---------------------------------------------------------------------- */

DumpVTK::DumpVTK(SPPARKS *spk, int narg, char **arg) : 
  DumpText(spk, narg, arg)
{
  // error check for allowed dump vtk output

  if (!multifile)
    error->all(FLERR,"Dump vtk must write one file per snapshot");
  if (multiproc) 
    error->all(FLERR,"Dump vtk cannot write multiple files per snapshot");
  if (binary) error->all(FLERR,"Dump vtk cannot write to binary files");

  // require a single per-site value

  if (size_one != 1) error->all(FLERR,"Dump vtk requires a single field");
  if (fields[0] == SITE || fields[0] == IARRAY) type = 0;
  else if (fields[0] == DARRAY) type = 1;
  else error->all(FLERR,"Dump vtk requires a field for a per-site value");

  // check for valid lattice

  if (!domain->lattice) 
    error->all(FLERR,"Dump vtk requires a lattice be defined");

  int latstyle = domain->lattice->style;
  int valid = 0;
  if (latstyle == LINE_2N || latstyle == SQ_4N || latstyle == SQ_8N || 
      latstyle == SC_6N || latstyle == SC_26N) valid = 1;
  if (!valid) 
    error->all(FLERR,"Dump vtk requires a regular lattice of sites");

  vtkflag = 0;
}

/* ---------------------------------------------------------------------- */

void DumpVTK::init_style()
{
  if (!vtkflag) error->all(FLERR,"Dump vtk requires dump modify vtk settings");
  if (!sort_flag || sortcol != 0) 
    error->all(FLERR,
               "Dump vtk requires ID-sorted output via dump modify command");
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_header(bigint ndump, double time)
{
  if (ndump != nx*ny*nz)
    error->all(FLERR,"Dump vtk output does not match "
               "dump modify vtk settings");

  if (me) return;

  fprintf(fp,"<?xml version=\"1.0\"\?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n");
  fprintf(fp,"<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" "
          "Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",nx,ny,nz);
  fprintf(fp,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",nx,ny,nz);
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData Scalars=\"Spin\">\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"Spin\" "
          "format=\"ascii\" RangeMin=\"%d\" RangeMax=\"%d\">\n",minval,maxval);
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_data(int n, double *buf)
{
  if (type == 0) {
    int m = 0;
    for (int i = 0; i < n; i++) {
      fprintf(fp,"%d ",static_cast<int> (buf[m]));
      m++;
    }
  } else {
    int m = 0;
    for (int i = 0; i < n; i++) {
      fprintf(fp,"%g ",buf[m]);
      m++;
    }
  }

  //fprintf(fp,"\n");
}

/* ---------------------------------------------------------------------- */

void DumpVTK::write_footer()
{
  if (me) return;

  fprintf(fp,"\n");
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
}

/* ---------------------------------------------------------------------- */

int DumpVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"vtk") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal dump_modify command");
    vtkflag = 1;
    nx = atoi(arg[1]);
    ny = atoi(arg[2]);
    nz = atoi(arg[3]);
    minval = atoi(arg[4]);
    maxval = atoi(arg[5]);
    return 6;
  }

  return 0;
}
