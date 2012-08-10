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

#include "spktype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dump_text.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

// customize by adding keyword to 1st enum

enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // also in dump_image
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,TAGINT};           // also in dump_image

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

DumpText::DumpText(SPPARKS *spk, int narg, char **arg) : Dump(spk, narg, arg)
{
  // allow for default "id site x y z" output if no args specified

  int def = 0;
  char **argcopy;
  
  if (narg == 4) {
    def = 1;
    narg = 9;
    argcopy = new char*[narg];
    argcopy[4] = (char *) "id";
    argcopy[5] = new char[6];
    strcpy(argcopy[5],"site");
    argcopy[6] = (char *) "x";
    argcopy[7] = (char *) "y";
    argcopy[8] = (char *) "z";
  } else argcopy = arg;

  // size_one may be shrunk below if additional optional args exist

  size_one = narg - 4;
  vtype = new int[size_one];
  vindex = new int[size_one];
  vformat = new char*[size_one];
  pack_choice = new FnPtrPack[size_one];

  // process attributes
  // ioptional = start of additional optional args
  // only dump image style processes optional args

  ioptional = parse_fields(narg,argcopy);

  if (ioptional < narg && strcmp(style,"image") != 0)
    error->all(FLERR,"Invalid attribute in dump text command");
  size_one = ioptional - 4;

  // setup vformat strings, one per field

  for (int i = 0; i < size_one; i++) {
    char *format;
    if (vtype[i] == INT) format = (char *) "%d ";
    else if (vtype[i] == DOUBLE) format = (char *) "%g ";
    else if (vtype[i] == TAGINT) format = (char *) TAGINT_FORMAT " ";
    int n = strlen(format) + 1;
    vformat[i] = new char[n];
    strcpy(vformat[i],format);
  }

  // setup column string
  // change "site" to "type" to be LAMMPS compatible

  int n = 0;
  for (int iarg = 4; iarg < narg; iarg++) n += strlen(argcopy[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 4; iarg < narg; iarg++) {
    if (strstr(argcopy[iarg],"site")) strcat(columns,"type");
    else strcat(columns,argcopy[iarg]);
    strcat(columns," ");
  }

  // delete argcopy if default output created

  if (def) {
    delete [] argcopy[5];
    delete [] argcopy;
  }

  // dump params

  iregion = -1;
  idregion = NULL;

  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;
  thresh_index = NULL;

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

  // setup function ptrs

  if (binary) header_choice = &DumpText::header_binary;
  else header_choice = &DumpText::header_text;

  if (binary) write_choice = &DumpText::write_binary;
  else write_choice = &DumpText::write_text;
}

/* ---------------------------------------------------------------------- */

DumpText::~DumpText()
{
  delete [] columns;

  delete [] idregion;
  memory->sfree(thresh_array);
  memory->sfree(thresh_op);
  memory->sfree(thresh_value);
  memory->sfree(thresh_index);

  memory->sfree(choose);
  memory->sfree(dchoose);
  memory->sfree(clist);

  delete [] vtype;
  delete [] vindex;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;
}

/* ---------------------------------------------------------------------- */

void DumpText::init_style()
{
  // error if propensity is dumped or used as threshold but doesn't exist
  // can't check until now, b/c input script may not have defined solver

  int flag = 0;
  for (int i = 0; i < size_one; i++)
    if (pack_choice[i] == &DumpText::pack_propensity) flag = 1;
  for (int i = 0; i < nthresh; i++)
    if (thresh_array[i] == PROPENSITY) flag = 1;
  if (flag && !solve)
    error->all(FLERR,"Dump requires propensity but no KMC solve performed");

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1) error->all(FLERR,"Region ID for dump text does not exist");
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpText::count()
{
  int i;

  // grow choose arrays if needed

  int nlocal = app->nlocal;
  if (nlocal > maxlocal) {
    maxlocal = nlocal;

    memory->sfree(choose);
    memory->sfree(dchoose);
    memory->sfree(clist);
    choose = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:choose");
    dchoose = (double *) 
      memory->smalloc(maxlocal*sizeof(double),"dump:dchoose");
    clist = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:clist");
  }

  // choose all local sites for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in region

  if (iregion >= 0) {
    Region *region = domain->regions[iregion];
    double **xyz = app->xyz;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(xyz[i][0],xyz[i][1],xyz[i][2]) == 0) 
	choose[i] = 0;
  }

  // un-choose if any threshhold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      if (thresh_array[ithresh] == ID) {
	tagint *id = app->id;
	for (i = 0; i < nlocal; i++) dchoose[i] = id[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == SITE) {
	int *site = app->iarray[0];
	for (i = 0; i < nlocal; i++) dchoose[i] = site[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == X) {
	ptr = &app->xyz[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
	ptr = &app->xyz[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
	ptr = &app->xyz[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == ENERGY) {
	if (latticeflag)
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = applattice->site_energy(i);
	else
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = appoff->site_energy(i);
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == PROPENSITY) {
	if (latticeflag)
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = applattice->site_propensity(i);
	else
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = appoff->site_propensity(i);
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == IARRAY) {
	int index = thresh_index[ithresh];
	if (latticeflag)
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = app->iarray[index][i];
	else
	  for (i = 0; i < nlocal; i++)
	    dchoose[i] = app->iarray[index][i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == DARRAY) {
	ptr = app->darray[thresh_index[ithresh]];
	nstride = 1;
      }
      
      // unselect sites that don't meet threshold criterion
      
      value = thresh_value[ithresh];
      
      if (thresh_op[ithresh] == LT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr >= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == LE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr > value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr <= value) choose[i] = 0;
      } else if (thresh_op[ithresh] == GE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr < value) choose[i] = 0;
      } else if (thresh_op[ithresh] == EQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr != value) choose[i] = 0;
      } else if (thresh_op[ithresh] == NEQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr == value) choose[i] = 0;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom
  
  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;
  
  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpText::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_header(int ndump, double time)
{
  if (multiproc) (this->*header_choice)(ndump,time);
  else if (me == 0) (this->*header_choice)(ndump,time);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_data(int n, double *buf)
{
  (this->*write_choice)(n,buf);
}

/* ---------------------------------------------------------------------- */

void DumpText::header_binary(int ndump, double time)
{
  fwrite(&idump,sizeof(int),1,fp);
  fwrite(&time,sizeof(double),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) {
    int one = 1;
    fwrite(&one,sizeof(int),1,fp);
  } else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpText::header_text(int ndump, double time)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d %10g\n",idump,time);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_binary(int n, double *buf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(buf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpText::write_text(int n, double *buf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT)
	fprintf(fp,vformat[j],static_cast<int> (buf[m]));
      else if (vtype[j] == DOUBLE)
	fprintf(fp,vformat[j],buf[m]);
      else if (vtype[j] == TAGINT) 
	fprintf(fp,vformat[j],static_cast<tagint> (buf[m]));
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpText::parse_fields(int narg, char **arg)
{
  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg-4;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpText::pack_id;
      vtype[i] = TAGINT;
    } else if (strcmp(arg[iarg],"site") == 0) {
      pack_choice[i] = &DumpText::pack_site;
      vtype[i] = INT;
      if (app->iarray == NULL)
	error->all(FLERR,"Dumping a quantity application does not support");
    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpText::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpText::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpText::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      pack_choice[i] = &DumpText::pack_energy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"propensity") == 0) {
      pack_choice[i] = &DumpText::pack_propensity;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (arg[iarg][0] == 'i') {
      pack_choice[i] = &DumpText::pack_iarray;
      vtype[i] = INT;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ninteger))
	error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;
    } else if (arg[iarg][0] == 'd') {
      pack_choice[i] = &DumpText::pack_darray;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&arg[iarg][1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ndouble))
	error->all(FLERR,"Invalid keyword in dump command");
      vindex[i]--;

    } else return iarg;
  }

  return narg;
}

/* ---------------------------------------------------------------------- */

int DumpText::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1) error->all(FLERR,"Dump_modify region ID does not exist");
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;

  } else if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
	memory->sfree(thresh_array);
	memory->sfree(thresh_op);
	memory->sfree(thresh_value);
	memory->sfree(thresh_index);
	thresh_array = NULL;
	thresh_op = NULL;
	thresh_value = NULL;
	thresh_index = NULL;
      }
      nthresh = 0;
      return 2;
    }
    
    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");
    
    // grow threshold arrays
    
    thresh_array = (int *)
      memory->srealloc(thresh_array,(nthresh+1)*sizeof(int),
		       "dump:thresh_array");
    thresh_op = (int *)
      memory->srealloc(thresh_op,(nthresh+1)*sizeof(int),
		       "dump:thresh_op");
    thresh_value = (double *)
      memory->srealloc(thresh_value,(nthresh+1)*sizeof(double),
		       "dump:thresh_value");
    thresh_index = (int *)
      memory->srealloc(thresh_index,(nthresh+1)*sizeof(int),
		       "dump:thresh_index");
    
    // set keyword type of threshold
    // customize by adding to if statement
    
    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    
    else if (strcmp(arg[1],"site") == 0) {
      if (app->iarray == NULL)
	error->all(FLERR,"Threshold for a quantity application does not support");
      thresh_array[nthresh] = SITE;
    }
    
    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;
    else if (strcmp(arg[1],"energy") == 0) thresh_array[nthresh] = ENERGY;
    else if (strcmp(arg[1],"propensity") == 0)
      thresh_array[nthresh] = PROPENSITY;
    
    // integer value = iN
    // double value = dN
    
    else if (arg[1][0] == 'i') {
      thresh_array[nthresh] = IARRAY;
      thresh_index[nthresh] = atoi(&arg[1][1]);
      if (thresh_index[nthresh] < 1 || 
	  thresh_index[nthresh] > app->ninteger)
	error->all(FLERR,"Threshold for a quantity application does not support");
      thresh_index[nthresh]--;
    } else if (arg[1][0] == 'd') {
      thresh_array[nthresh] = DARRAY;
      thresh_index[nthresh] = atoi(&arg[1][1]);
      if (thresh_index[nthresh] < 1 || 
	  thresh_index[nthresh] > app->ndouble)
	error->all(FLERR,"Threshold for a quantity application does not support");
      thresh_index[nthresh]--;
      
    } else error->all(FLERR,"Invalid dump_modify threshold operator");
    
    // set operation type of threshold
    
    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all(FLERR,"Invalid dump_modify threshold operator");
    
    // set threshold value
    
    thresh_value[nthresh] = atof(arg[3]);
    
    nthresh++;
    return 4;
  }

  return 0;
}


/* ----------------------------------------------------------------------
   one method for every keyword dump can output
   the site quantity is packed into buf starting at n with stride size_one
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void DumpText::pack_id(int n)
{
  tagint *id = app->id;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = id[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_site(int n)
{
  int *site = app->iarray[0];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = site[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_x(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = xyz[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_y(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = xyz[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_z(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = xyz[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_energy(int n)
{
  if (latticeflag) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = applattice->site_energy(clist[i]);
      n += size_one;
    }
  } else {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = appoff->site_energy(clist[i]);
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_propensity(int n)
{
  if (latticeflag) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = applattice->site_propensity(clist[i]);
      n += size_one;
    }
  } else {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = appoff->site_propensity(clist[i]);
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_iarray(int n)
{
  int *ivec = app->iarray[vindex[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = ivec[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpText::pack_darray(int n)
{
  double *dvec = app->darray[vindex[n]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = dvec[clist[i]];
    n += size_one;
  }
}
