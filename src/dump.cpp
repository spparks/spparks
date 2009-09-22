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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "dump.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

// customize by adding keyword to 1st enum

enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE};

#define MAXLINE 1024
#define DEFAULT "id site x y z"

/* ---------------------------------------------------------------------- */

Dump::Dump(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // parse dump args

  if (narg < 3) error->all("Illegal dump command");

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  delta = atof(arg[1]);
  if (delta <= 0.0) error->all("Illegal dump command");

  n = strlen(arg[2]) + 1;
  filename = new char[n];
  strcpy(filename,arg[2]);

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  compressed = 0;
  binary = 0;
  multifile = 0;
  multiproc = 0;

  char *ptr;
  if (ptr = strchr(filename,'%')) {
    multiproc = 1;
    char *extend = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(extend,"%s%d%s",filename,me,ptr+1);
    delete [] filename;
    n = strlen(extend) + 1;
    filename = new char[n];
    strcpy(filename,extend);
    delete [] extend;
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  } else
    error->all("Dump command can only be used for spatial applications");

  // parse fields
  // use DEFAULT if fields not listed

  columns = new char[MAXLINE];

  if (narg-3 == 0) {
    strcpy(columns,DEFAULT);
    size_one = 5;
  } else {
    columns[0] = '\0';
    for (int iarg = 3; iarg < narg; iarg++) {
      strcat(columns,arg[iarg]);
      strcat(columns," ");
    }
    columns[strlen(columns)-1] = '\0';
    size_one = narg-3;
  }

  // parse columns string into keywords
  // make a copy so columns can be output in dump header

  vtype = new int[size_one];
  vindex = new int[size_one];
  vformat = new char*[size_one];
  pack_choice = new FnPtrPack[size_one];

  n = strlen(columns) + 1;
  char *line = new char[n];
  strcpy(line,columns);

  int i = 0;
  char *word = strtok(line," \0");
  while (word) {
    if (strcmp(word,"id") == 0) {
      pack_choice[i] = &Dump::pack_id;
      vtype[i] = INT;
    } else if (strcmp(word,"site") == 0) {
      pack_choice[i] = &Dump::pack_site;
      vtype[i] = INT;
      if (app->iarray == NULL)
	error->all("Dumping a quantity application does not support");
    } else if (strcmp(word,"x") == 0) {
      pack_choice[i] = &Dump::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"y") == 0) {
      pack_choice[i] = &Dump::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"z") == 0) {
      pack_choice[i] = &Dump::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"energy") == 0) {
      pack_choice[i] = &Dump::pack_energy;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"propensity") == 0) {
      pack_choice[i] = &Dump::pack_propensity;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (word[0] == 'i') {
      pack_choice[i] = &Dump::pack_iarray;
      vtype[i] = INT;
      vindex[i] = atoi(&word[1]);
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ninteger))
	error->all("Invalid keyword in dump command");
      vindex[i]--;
    } else if (word[0] == 'd') {
      pack_choice[i] = &Dump::pack_darray;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&word[1]) - 1;
      if (latticeflag && (vindex[i] < 1 || vindex[i] > app->ndouble))
	error->all("Invalid keyword in dump command");
      vindex[i]--;

    } else error->all("Invalid keyword in dump command");

    word = strtok(NULL," \0");
    i++;
  }

  delete [] line;

  // setup vformat strings, one per field

  for (int i = 0; i < size_one; i++) {
    char *format;
    if (vtype[i] == INT) format = "%d ";
    else if (vtype[i] == DOUBLE) format = "%g ";
    int n = strlen(format) + 1;
    vformat[i] = new char[n];
    strcpy(vformat[i],format);
  }

  // in columns, change "site" to "type" to be LAMMPS compatible

  if (strstr(columns,"site")) {
    char *ptr1 = strstr(columns,"site");
    char *ptr2 = ptr1 + strlen("site");
    *ptr1 = '\0';
    sprintf(columns,"%s%s%s",columns,"type",ptr2);
  }

  // dump params

  if (latticeflag) nlocal = applattice->nlocal;

  boxxlo = app->boxxlo;
  boxxhi = app->boxxhi;
  boxylo = app->boxylo;
  boxyhi = app->boxyhi;
  boxzlo = app->boxzlo;
  boxzhi = app->boxzhi;

  flush_flag = 1;
  logfreq = 0;
  delay = 0.0;

  idump = 0;

  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;
  thresh_index = NULL;

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;

  maxbuf = 0;
  buf = NULL;

  // one-time file open

  if (multifile == 0) openfile();

  // setup function ptrs

  if (binary) header_choice = &Dump::header_binary;
  else header_choice = &Dump::header_text;

  if (binary) write_choice = &Dump::write_binary;
  else write_choice = &Dump::write_text;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] filename;
  delete [] columns;

  memory->sfree(thresh_array);
  memory->sfree(thresh_op);
  memory->sfree(thresh_value);
  memory->sfree(thresh_index);

  memory->sfree(choose);
  memory->sfree(dchoose);
  memory->sfree(buf);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }

  delete [] vtype;
  delete [] vindex;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  // error if propensity is dumped or used as threshold but doesn't exist
  // can't check until now, b/c input script may not have defined solver

  int flag = 0;
  for (int i = 0; i < size_one; i++)
    if (pack_choice[i] == &Dump::pack_propensity) flag = 1;
  for (int i = 0; i < nthresh; i++)
    if (thresh_array[i] == PROPENSITY) flag = 1;
  if (flag && !solve)
    error->all("Dump requires propensity but no KMC solve performed");
}

/* ---------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delta") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      delta = atof(arg[iarg+1]);
      if (delta <= 0.0) error->all("Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all("Illegal dump_modify command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (nrepeat < 0) error->all("Illegal dump_modify command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      delay = atof(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"thresh") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"none") == 0) {
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
	iarg += 2;
	continue;
      }
      
      if (iarg+4 > narg) error->all("Illegal dump_modify command");
      
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
      
      if (strcmp(arg[iarg+1],"id") == 0) thresh_array[nthresh] = ID;

      else if (strcmp(arg[iarg+1],"site") == 0) {
	if (app->iarray == NULL)
	  error->all("Threshold for a quantity application does not support");
	thresh_array[nthresh] = SITE;
      }

      else if (strcmp(arg[iarg+1],"x") == 0) thresh_array[nthresh] = X;
      else if (strcmp(arg[iarg+1],"y") == 0) thresh_array[nthresh] = Y;
      else if (strcmp(arg[iarg+1],"z") == 0) thresh_array[nthresh] = Z;

      else if (strcmp(arg[iarg+1],"energy") == 0)
	thresh_array[nthresh] = ENERGY;
      else if (strcmp(arg[iarg+1],"propensity") == 0)
	thresh_array[nthresh] = PROPENSITY;
      
      // integer value = iN
      // double value = dN

      else if (arg[iarg+1][0] == 'i') {
	thresh_array[nthresh] = IARRAY;
	thresh_index[nthresh] = atoi(&arg[iarg+1][1]);
	if (thresh_index[nthresh] < 1 || 
	    thresh_index[nthresh] > app->ninteger)
	  error->all("Threshold for a quantity application does not support");
	thresh_index[nthresh]--;
      } else if (arg[iarg+1][0] == 'd') {
	thresh_array[nthresh] = DARRAY;
	thresh_index[nthresh] = atoi(&arg[iarg+1][1]);
	if (thresh_index[nthresh] < 1 || 
	    thresh_index[nthresh] > app->ndouble)
	  error->all("Threshold for a quantity application does not support");
	thresh_index[nthresh]--;

      } else error->all("Invalid dump_modify threshold operator");
      
      // set operation type of threshold
      
      if (strcmp(arg[iarg+2],"<") == 0) thresh_op[nthresh] = LT;
      else if (strcmp(arg[iarg+2],"<=") == 0) thresh_op[nthresh] = LE;
      else if (strcmp(arg[iarg+2],">") == 0) thresh_op[nthresh] = GT;
      else if (strcmp(arg[iarg+2],">=") == 0) thresh_op[nthresh] = GE;
      else if (strcmp(arg[iarg+2],"==") == 0) thresh_op[nthresh] = EQ;
      else if (strcmp(arg[iarg+2],"!=") == 0) thresh_op[nthresh] = NEQ;
      else error->all("Invalid dump_modify threshold operator");
      
      // set threshold value
      
      thresh_value[nthresh] = atof(arg[iarg+3]);
      nthresh++;
      iarg += 4;

    } else error->all("Illegal dump_modify command");
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of site values as atom coords
------------------------------------------------------------------------- */

void Dump::write(double time)
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // nmine = # of dump lines this proc will contribute to dump
  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int nmine = count();

  int ntotal,nmax;
  if (multiproc) nmax = nmine;
  else {
    MPI_Allreduce(&nmine,&ntotal,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nmine,&nmax,1,MPI_INT,MPI_MAX,world);
  }

  // write timestep header

  if (multiproc) write_header(nmine,time);
  else write_header(ntotal,time);

  idump++;

  // grow communication buffer if necessary

  if (nmax*size_one > maxbuf) {
    maxbuf = nmax*size_one;
    memory->sfree(buf);
    buf = (double *) memory->smalloc(maxbuf*sizeof(double),"dump:buf");
  }

  // pack my data into buf

  pack();

  // multiproc = 1 = each proc writes own data to own file 
  // multiproc = 0 = all procs write to one file thru proc 0
  //   proc 0 pings each proc, receives it's data, writes to file
  //   all other procs wait for ping, send their data to proc 0

  if (multiproc) write_data(nmine,buf);
  else {
    int tmp,nlines;
    MPI_Status status;
    MPI_Request request;

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
	if (iproc) {
	  MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	  nlines /= size_one;
	} else nlines = nmine;

	write_data(nlines,buf);
      }
      if (flush_flag) fflush(fp);
      
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,nmine*size_one,MPI_DOUBLE,0,0,world);
    }
  }

  // if file per timestep, close file

  if (multifile) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

int Dump::count()
{
  int i;

  // update nlocal if off-lattice app

  if (latticeflag == 0) nlocal = appoff->nlocal;

  // insure choose arrays are big enough
  // initialize choose[] to select all sites if reallocate

  if (nlocal > maxlocal) {
    maxlocal = nlocal;

    memory->sfree(choose);
    choose = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:choose");
    for (int i = 0; i < nlocal; i++) choose[i] = 1;

    if (nthresh) {
      memory->sfree(dchoose);
      dchoose = (double *) 
	memory->smalloc(maxlocal*sizeof(double),"dump:dchoose");
    }
  }

  // no thresholds, select all sites

  if (nthresh == 0) return nlocal;

  // unselect a site if any threshold criterion isn't met
  // dchoose stores the site value to compare threshold value to

  double *ptr;
  double value;
  int nstride;

  int nmine = nlocal;
    
  for (int ithresh = 0; ithresh < nthresh; ithresh++) {

    // customize by adding to if statement

    if (thresh_array[ithresh] == ID) {
      int *id = app->id;
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
	if (choose[i] && *ptr >= value) {
	  choose[i] = 0;
	  nmine--;
	}
    } else if (thresh_op[ithresh] == LE) {
      for (i = 0; i < nlocal; i++, ptr += nstride)
	if (choose[i] && *ptr > value) {
	  choose[i] = 0;
	  nmine--;
	}
    } else if (thresh_op[ithresh] == GT) {
      for (i = 0; i < nlocal; i++, ptr += nstride)
	if (choose[i] && *ptr <= value) {
	  choose[i] = 0;
	  nmine--;
	}
    } else if (thresh_op[ithresh] == GE) {
      for (i = 0; i < nlocal; i++, ptr += nstride)
	if (choose[i] && *ptr < value) {
	  choose[i] = 0;
	  nmine--;
	}
    } else if (thresh_op[ithresh] == EQ) {
      for (i = 0; i < nlocal; i++, ptr += nstride)
	if (choose[i] && *ptr != value) {
	  choose[i] = 0;
	  nmine--;
	}
    } else if (thresh_op[ithresh] == NEQ) {
      for (i = 0; i < nlocal; i++, ptr += nstride)
	if (choose[i] && *ptr == value) {
	  choose[i] = 0;
	  nmine--;
	}
    }
  }

  return nmine;
}

/* ---------------------------------------------------------------------- */

void Dump::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void Dump::write_header(int ndump, double time)
{
  if (multiproc) (this->*header_choice)(ndump,time);
  else if (me == 0) (this->*header_choice)(ndump,time);
}

/* ---------------------------------------------------------------------- */

void Dump::write_data(int n, double *buf)
{
  (this->*write_choice)(n,buf);
}

/* ---------------------------------------------------------------------- */

void Dump::header_binary(int ndump, double time)
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

void Dump::header_text(int ndump, double time)
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

void Dump::write_binary(int n, double *buf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(buf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void Dump::write_text(int n, double *buf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (buf[m]));
      else fprintf(fp,vformat[j],buf[m]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // if one file per timestep, replace '*' with current timestep

  char *filecurrent;
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    sprintf(filecurrent,"%s%d%s",filename,idump,ptr+1);
    *ptr = '*';
  }

  // open one file on proc 0 or file on every proc

  if (me == 0 || multiproc) {
    if (compressed) {
#ifdef SPPARKS_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
      fp = popen(gzip,"w");
#else
      error->one("Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == NULL) error->one("Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ----------------------------------------------------------------------
   one method for every keyword dump can output
   the site quantity is packed into buf starting at n with stride size_one
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Dump::pack_id(int n)
{
  int *id = app->id;

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = id[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_site(int n)
{
  int *site = app->iarray[0];

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = site[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_x(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = xyz[i][0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_y(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = xyz[i][1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_z(int n)
{
  double **xyz = app->xyz;

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = xyz[i][2];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_energy(int n)
{
  if (latticeflag) {
    for (int i = 0; i < nlocal; i++) {
      if (choose[i]) {
	buf[n] = applattice->site_energy(i);
	n += size_one;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (choose[i]) {
	buf[n] = appoff->site_energy(i);
	n += size_one;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_propensity(int n)
{
  if (latticeflag) {
    for (int i = 0; i < nlocal; i++) {
      if (choose[i]) {
	buf[n] = applattice->site_propensity(i);
	n += size_one;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (choose[i]) {
	buf[n] = appoff->site_propensity(i);
	n += size_one;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_iarray(int n)
{
  int *ivec = app->iarray[vindex[n]];

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = ivec[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::pack_darray(int n)
{
  double *dvec = app->darray[vindex[n]];

  for (int i = 0; i < nlocal; i++) {
    if (choose[i]) {
      buf[n] = dvec[i];
      n += size_one;
    }
  }
}
