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

#include <cstdio>
#include "spktype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "unistd.h"
#include "dump_stitch.h"
#include "app.h"
#include "app_lattice.h"
#include "domain.h"
#include "error.h"
#include <sstream>
#include <iomanip>
#include <iostream>

#include "stitch.h"

using namespace SPPARKS_NS;

// customize by adding keyword to 1st enum

enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // in other dump files
enum{INT,DOUBLE,TAGINT};           // in other dump files

/* ---------------------------------------------------------------------- */

DumpStitch::DumpStitch(SPPARKS *spk, int narg, char **arg) : Dump(spk, narg, arg)
{
  if (multifile || multiproc)
    error->all(FLERR,"Dump stitch cannot use '*' or '%' wildcards");

  if (narg < 4) error->all(FLERR,"Illegal dump stitch command");

  size_one = narg - 4;
  fields = new int[size_one];
  vtype = new int[size_one];
  vindex = new int[size_one];
  stitch_field_ids = new int64_t[size_one];

  parse_fields(narg,arg);
  create_stitch_field_ids();

  // error check than only integer fields were specified

  for (int i = 0; i < size_one; i++)
    if (INT!=vtype[i] && DOUBLE != vtype[i])
      error->all(FLERR,"Dump stitch only supports 'integer' or 'double'  values.");
}

/* ---------------------------------------------------------------------- */

DumpStitch::~DumpStitch()
{
  delete [] fields;
  delete [] vtype;
  delete [] vindex;
  delete [] stitch_field_ids;

  int err = stitch_close(&stitch_file);
}

/* ---------------------------------------------------------------------- */

void DumpStitch::init_style()
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Dump stitch only allowed for on-lattice apps");

  applattice = (AppLattice*) app;

  if (!applattice->simple) 
    error->all(FLERR,"Dump stitch only allowed for simple lattices");
  
   // error if propensity is dumped but doesn't exist
   // can't check until now, b/c input script may not have defined solver

   int flag = 0;
   for (int i = 0; i < size_one; i++)
     if (fields[i] == PROPENSITY) flag = 1;
   if (flag && !solve)
     error->all(FLERR,"Dump requires proensity but no KMC solve performed");
}

/* ---------------------------------------------------------------------- */

void DumpStitch::create_stitch_field_ids(){

   // open Stitch database file
   int err = stitch_open(&stitch_file,MPI_COMM_WORLD,filename);
   // TODO: process err
   //

   int64_t field_id=-1;
   for (int i = 0; i < size_one; i++) {
      char label[8];
      if (SITE==fields[i]) {
         err = stitch_query_field(stitch_file,"site",&field_id);
         // TODO: process err
         //
         if(STITCH_FIELD_NOT_FOUND_ID==field_id){
             // Specify value type of value and its undefined value
             union StitchTypesUnion v;v.i32=-1;
             // Scalar value
             int32_t scalar=1;
             err = stitch_create_field(stitch_file,"site",STITCH_INT32,
				       v,scalar,&field_id);
             // TODO: process err
             //
         }
      }
      else if (IARRAY==fields[i]) {
         char label[8];
         sprintf(label,"i%d",vindex[i]+1);
         err = stitch_query_field(stitch_file,label,&field_id);
         // TODO: process err
         //
         if(STITCH_FIELD_NOT_FOUND_ID==field_id){
             // Specify value type of value and its undefined value
             union StitchTypesUnion v;v.i32=-1;
             // Scalar value
             int32_t scalar=1;
             err = stitch_create_field(stitch_file,label,STITCH_INT32,
				       v,scalar,&field_id);
             // TODO: process err
             //
         }
      }
      else if (DARRAY==fields[i] || X==fields[i] || Y==fields[i] || Z==fields[i]) {
         char label[8];
         sprintf(label,"d%d",vindex[i]+1);
         err = stitch_query_field(stitch_file,label,&field_id);
         // TODO: process err
         //
         if(STITCH_FIELD_NOT_FOUND_ID==field_id){
             // Specify value type of value and its undefined value
             union StitchTypesUnion v;v.f64=-1;
             // Scalar value
             int32_t scalar=1;
             err = stitch_create_field (stitch_file,label,STITCH_FLOAT64,
					v,scalar,&field_id);
             // TODO: process err
             //
         }
      } else {
         error->all(FLERR,"DumpStitch::create_stitch_field_ids(); "
		    "Invalid specified field.");
      }
      stitch_field_ids[i]=field_id;
   }
}

/* ----------------------------------------------------------------------
   dump a snapshot of site values
------------------------------------------------------------------------- */

void DumpStitch::write(double time)
{
  int block[6];
  block[0] = applattice->xlo_me_simple;
  block[1] = applattice->xhi_me_simple+1;
  block[2] = applattice->ylo_me_simple;
  block[3] = applattice->yhi_me_simple+1;
  block[4] = applattice->zlo_me_simple;
  block[5] = applattice->zhi_me_simple+1;
  int xlo = block[0];
  int xhi = block[1];
  int ylo = block[2];
  int yhi = block[3];
  int zlo = block[4];
  int zhi = block[5];

#ifdef LOG_STITCH
  {
     int my_rank, num_procs;
     MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
     std::size_t num_digits=1;
     if(num_procs>=10)
       num_digits=2;
     else if(num_procs>=100)
       num_digits=3;

     std::ostringstream fname;
     fname << "write.blocks.np" << std::setw(num_digits) << std::setfill('0') << num_procs << ".";
     // fills rank leading '0' 
     fname << std::setw(num_digits) << std::setfill('0') << my_rank << ".dat";
     FILE* fp = std::fopen(fname.str().c_str(), "a");
     if(!fp){
        error->all(FLERR,"dump_stitch.cpp; in function 'write' file opening failed.");
     }
     const char* fmt="%10.4f,%4d,%4d,%4d,%4d,%4d,%4d\n";
     fprintf(fp,fmt,time,xlo,xhi,ylo,yhi,zlo,zhi);
     std::fclose(fp);
  }
#endif


  int *idata = NULL;
  double *real_data = NULL;
  int err,write_flag = -1;
  int64_t field_id=-1;

  for (int i = 0; i < size_one; i++) {
     field_id=stitch_field_ids[i];
     if (SITE==fields[i]) {
       idata = app->iarray[0];
       //For write:
       //write_flag == True when new time step created
       //write_flag == False when new time step NOT created
       //int stitch_write_block_int32 (const StitchFile * file, int64_t field_id, double * time, int32_t * bb, const int32_t * buffer, int32_t * write_flag);
       err = stitch_write_block_int32(stitch_file,field_id,&time,
				      block,idata,&write_flag);
       //printf("1:dump_stitch; time=%3.1f, write_flag=%d\n",time,write_flag);
       //printf("2: block xlo, xhi, ylo, yhi, zlo, zli = %5d,%5d,%5d,%5d,%5d,%5d\n",xlo,xhi,ylo,yhi,zlo,zhi);
     }
     else if (IARRAY==fields[i]) {
       //std::cout << "DumpStitch::write(time=" << time << "); IARRAY field_id="<<field_id<< std::endl;
       idata = app->iarray[vindex[i]];
       err = stitch_write_block_int32(stitch_file,field_id,&time,
				      block,idata,&write_flag);
     }
     else if (DARRAY==fields[i]) {
       //std::cout << "DumpStitch::write(time=" << time << "); DARRAY field_id="<<field_id<< std::endl;
       real_data = app->darray[vindex[i]];
       err = stitch_write_block_float64(stitch_file,field_id,&time,
					block,real_data,&write_flag);
     }
     // TODO: process err
     //

    //int rank=-1;
    //err=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    //printf("START dump stitch rank = %d\n",rank);
    //MPI_Barrier(MPI_COMM_WORLD);
    //err = stitch_write_block(stitch_file,&time,block,idata,&write_flag);
    //printf("1:dump_stitch; time=%3.1f, write_flag=%d\n",time,write_flag);
    //printf("2: block xlo, xhi, ylo, yhi, zlo, zli = %5d,%5d,%5d,%5d,%5d,%5d\n",xlo,xhi,ylo,yhi,zlo,zhi);
    //printf("FINISH dump stitch rank = %d\n",rank);
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("FINISH MPI_Barrier dump stitch rank = %d\n",rank);
  }
}

/* ---------------------------------------------------------------------- */

int DumpStitch::parse_fields(int narg, char **arg)
{
  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"id") == 0) {
      fields[i] = ID;
      vtype[i] = TAGINT;
    } else if (strcmp(arg[iarg],"site") == 0) {
      fields[i] = SITE;
      vtype[i] = INT;
      if (app->iarray == NULL)
         error->all(FLERR,"Dumping a quantity application does not support");
    } else if (strcmp(arg[iarg],"x") == 0) {
      fields[i] = X;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      fields[i] = Y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      fields[i] = Z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      fields[i] = ENERGY;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"propensity") == 0) {
      fields[i] = PROPENSITY;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (arg[iarg][0] == 'i') {
      fields[i] = IARRAY;
      vtype[i] = INT;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ninteger))
         error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;
    } else if (arg[iarg][0] == 'd') {
      fields[i] = DARRAY;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ndouble))
         error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;

    } else return iarg;
  }

  return narg;
}
