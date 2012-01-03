/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
 
   This class added by Eric Homer (copied after diag_energy.cpp).

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eric Homer (Sandia)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_array.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "comm_lattice.h"
#include "comm_off_lattice.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum {MEAN,SUM,MIN,MAX};

/* ---------------------------------------------------------------------- */

DiagArray::DiagArray(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass == App::LATTICE) latticeflag = 1;
  else if (app->appclass == App::OFF_LATTICE) latticeflag = 0;
  else error->all(FLERR,"Diag style incompatible with app style");
  
  // check the number of arguments

  if (((narg-1) % 2) != 0)
    error->all(FLERR,"Invalid diag_style command");
  
  // list the number of output variables

  nvals = (narg-1)/2;
  
  memory->create(index,nvals,"diag_array::index");
  memory->create(index_double,nvals,"diag_array::index_double");
  memory->create(vals,nvals,"diag_array::vals");
   memory->create(diag_method,nvals,"diag_array::diag_method");
  
  for (int i=0; i<nvals; i++) {

    // shift to zero based array
    
    index[i] = atoi(&arg[1+2*i][1])-1;
    
    // determine whether integer or double and check

    if (arg[1+2*i][0] == 'i') {
      index_double[i]=false;
      if (index[i]>= app->ninteger) 
        error->all(FLERR,"Invalid diag_style command");
    } else if (arg[1+2*i][0] == 'd') {
      index_double[i]=true;
      if (index[i]>= app->ndouble) 
        error->all(FLERR,"Invalid diag_style command");
    } else error->all(FLERR,"Invalid diag_style command");
    
    // determine which stat is requested

    if ( strcmp(arg[1+2*i+1],"mean") == 0)
      diag_method[i]=MEAN;
    else if ( strcmp(arg[1+2*i+1],"sum") == 0)
      diag_method[i]=SUM;
    else if ( strcmp(arg[1+2*i+1],"min") == 0)
      diag_method[i]=MIN;
    else if ( strcmp(arg[1+2*i+1],"max") == 0)
      diag_method[i]=MAX;
    else
      error->all(FLERR,"Invalid diag_style command");
  }
}

/* ---------------------------------------------------------------------- */

DiagArray::~DiagArray()
{
  memory->destroy(index);
  memory->destroy(index_double);
  memory->destroy(vals);
  memory->destroy(diag_method);
}

/* ---------------------------------------------------------------------- */

void DiagArray::init()
{
  if (latticeflag) applattice = (AppLattice *) app;
  else appofflattice = (AppOffLattice *) app;
}

/* ---------------------------------------------------------------------- */

void DiagArray::compute()
{
  int nlocal = app->nlocal;
  int nlocal_all;
  
  MPI_Allreduce(&nlocal,&nlocal_all,1,MPI_INT,MPI_SUM,world);
  
  double **dptr;
  int **iptr;
  
  double dtmp,dtmp_all;
  int itmp,itmp_all;
  int ii;
  
  
  if (latticeflag) {
    applattice->comm->all();
    iptr = applattice->iarray;
    dptr = applattice->darray;
  } else {
    appofflattice->comm->all();
    iptr = appofflattice->iarray;
    dptr = appofflattice->darray;
  }
  
  for (int i=0; i<nvals; i++) {
    
    ii=index[i];
    
    if (index_double[i]) {

      // calculate the diagnostic for double arrays
      
      if (diag_method[i]==MIN) {
        dtmp = dptr[ii][0];
        dtmp_all = dptr[ii][0];
        for (int j = 1; j < nlocal; j++) 
          if (dptr[ii][j] < dtmp)
            dtmp=dptr[ii][j];
        MPI_Allreduce(&dtmp,&dtmp_all,1,MPI_DOUBLE,MPI_MIN,world);
        vals[i]=dtmp_all;
      } else if (diag_method[i]==MAX) {
        dtmp = dptr[ii][0];
        dtmp_all = dptr[ii][0];
        for (int j = 1; j < nlocal; j++) 
          if (dptr[ii][j] > dtmp)
            dtmp=dptr[ii][j];
        MPI_Allreduce(&dtmp,&dtmp_all,1,MPI_DOUBLE,MPI_MAX,world);
        vals[i]=dtmp_all;
      } else {
        dtmp = 0.0;
        for (int j = 0; j < nlocal; j++) dtmp += dptr[ii][j];
        MPI_Allreduce(&dtmp,&dtmp_all,1,MPI_DOUBLE,MPI_SUM,world);
      
        if (diag_method[i]==MEAN)
          vals[i] = dtmp_all / double(nlocal_all);
        else if (diag_method[i]==SUM)
          vals[i] = dtmp_all;
      }
    } else {

      // calculate the diagnostic for integer arrays
      
      if (diag_method[i]==MIN) {
        itmp = iptr[ii][0];
        itmp_all = iptr[ii][0];
        for (int j = 1; j < nlocal; j++) 
          if (iptr[ii][j] < itmp)
            itmp=iptr[ii][j];
        MPI_Allreduce(&itmp,&itmp_all,1,MPI_INT,MPI_MIN,world);
        vals[i]= double(itmp_all);
      } else if (diag_method[i]==MAX) {
        itmp = iptr[ii][0];
        itmp_all = iptr[ii][0];
        for (int j = 1; j < nlocal; j++) 
          if (iptr[ii][j] > itmp)
            itmp=iptr[ii][j];
        MPI_Allreduce(&itmp,&itmp_all,1,MPI_INT,MPI_MAX,world);
        vals[i]= double(itmp_all);
      } else {
        itmp = 0;
        for (int j = 0; j < nlocal; j++) itmp += iptr[ii][j];
        MPI_Allreduce(&itmp,&itmp_all,1,MPI_INT,MPI_SUM,world);
        
        if (diag_method[i]==MEAN)
          vals[i] = double(itmp_all) / double(nlocal_all);
        else if (diag_method[i]==SUM)
          vals[i] = double(itmp_all);
      }  
    }
  } 
}

/* ---------------------------------------------------------------------- */

void DiagArray::stats(char *strtmp)
{
  for (int i=0; i<nvals; i++) {
    sprintf(strtmp," %10g",vals[i]);
    strtmp += strlen(strtmp);
  }
}

/* ---------------------------------------------------------------------- */

void DiagArray::stats_header(char *strtmp)
{
  char str1[10],str2[10];
  
  for (int i=0; i<nvals; i++) {
    
    if (diag_method[i]==MEAN)
      sprintf(str1,"%s","mean");
    else if (diag_method[i]==SUM)
      sprintf(str1,"%s","sum");
    else if (diag_method[i]==MIN)
      sprintf(str1,"%s","min");
    else if (diag_method[i]==MAX)
      sprintf(str1,"%s","max");
    
    // shift index back to 1-based array

    if (index_double[i])
      sprintf(str2,"%s(d%d)",str1,index[i]+1);  
    else 
      sprintf(str2,"%s(i%d)",str1,index[i]+1);
    
    sprintf(strtmp," %10s",str2);
    strtmp += strlen(strtmp);
  }
}
