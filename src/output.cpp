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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "app.h"
#include "app_lattice.h"
#include "diag.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{GENERAL,LATTICE};     // same as in app.h
enum{INT,DOUBLE};

#define MAXLINE 1024
#define DEFAULT "id lattice x y z"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Output::Output(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  stats_delta = 0.0;
  stats_logfreq = 0;
  stats_delay = 0.0;

  dump_delta = 0.0;
  dump_logfreq = 0;
  dump_delay = 0.0;
  idump = 0;

  ndiags = 0;
  diaglist = NULL;

  size_one = 0;
  fp = NULL;
  vtype = NULL;
  vindex = NULL;
  vformat = NULL;
  pack_choice = NULL;
  buf = NULL;
  mask = NULL;

  mask_flag = 0;
}

/* ---------------------------------------------------------------------- */

Output::~Output()
{
  if (me == 0 && fp) fclose(fp);

  for (int i = 0; i < ndiags; i++) delete diaglist[i];
  memory->sfree(diaglist);

  delete [] vtype;
  delete [] vindex;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;

  memory->sfree(mask);
}

/* ---------------------------------------------------------------------- */

void Output::init(double time)
{
  // error if dump is defined and propensity is output but doesn't exist

  if (dump_delta > 0.0) {
    int flag = 0;
    for (int i = 0; i < size_one; i++)
      if (pack_choice[i] == &Output::pack_propensity) flag = 1;
    if (flag && !solve)
      error->all("Dumping propensity but no KMC solve performed");
  }

  // initialize all diagnostics

  for (int i = 0; i < ndiags; i++) diaglist[i]->init();
}

/* ----------------------------------------------------------------------
   called before every run
   perform stats output
   set next output time for all kinds of output
   return tnext = next time any output is needed
------------------------------------------------------------------------- */

double Output::setup(double time)
{
  // initial dump file snapshot
  // needs to happen in setup() in case propensity is output
  // app not ready to compute propensities until setup_app() is called

  if (dump_delta > 0.0 && idump == 0) dump(time);

  // set next dump time

  if (dump_delta > 0.0) {
    dump_time = next_time(time,dump_logfreq,dump_delta,
			  dump_nrepeat,dump_scale,dump_delta);
  } else dump_time = app->stoptime;

  // if stats drives output, perform diagnostics
  // else set next diagnostic times

  double diag_time = app->stoptime;
  for (int i = 0; i < ndiags; i++) {
    if (diaglist[i]->stats_flag) diaglist[i]->compute();
    else {
      diaglist[i]->diag_time = 
	next_time(time,diaglist[i]->diag_logfreq,diaglist[i]->diag_delta,
		  diaglist[i]->diag_nrepeat,diaglist[i]->diag_scale,
		  diaglist[i]->diag_delta);
      diag_time = MIN(diag_time,diaglist[i]->diag_time);
    }
  }

  // perform stats output
  // set next stats time

  stats_header();
  stats(0);
  if (stats_delta > 0.0)
    stats_time = next_time(time,stats_logfreq,stats_delta,
			   stats_nrepeat,stats_scale,stats_delta);
  else stats_time = app->stoptime;

  double tnext = app->stoptime;
  tnext = MIN(tnext,dump_time);
  tnext = MIN(tnext,diag_time);
  tnext = MIN(tnext,stats_time);
  return tnext;
}

/* ---------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);
  if (stats_delta < 0.0) error->all("Illegal stats command");

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      if (stats_delta == 0.0) error->all("Illegal stats command");
      stats_logfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	stats_nrepeat = atoi(arg[iarg]);
	if (stats_nrepeat < 1) error->all("Illegal stats command");
	iarg++;
	stats_scale = atof(arg[iarg]);
	if (stats_nrepeat*stats_delta > stats_scale)
	  error->all("Illegal stats command");
      } else error->all("Illegal stats command");
    } else error->all("Illegal stats command");
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

void Output::set_dump(int narg, char **arg)
{
  // determine correct kind of app pointer

  if (app->appclass == LATTICE)
    applattice = (AppLattice *) app;
  else error->all("Cannot use dump with off-lattice app");

  // parse dump args

  if (narg < 2) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  if (me == 0) {
    if (fp) fclose(fp);
    fp = fopen(arg[1],"w");
    if (!fp) error->one("Cannot open dump file");
  }

  int iarg = 2;
  while (iarg < narg) {

  if (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      dump_logfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	dump_nrepeat = atoi(arg[iarg]);
	if (dump_nrepeat < 1) error->all("Illegal dump command");
	iarg++;
	dump_scale = atof(arg[iarg]);
	if (dump_nrepeat*dump_delta > dump_scale)
	  error->all("Illegal dump command");
	iarg++;
      } else {
	error->all("Illegal dump command");
      }
    } else if (strcmp(arg[iarg],"delay") == 0) {
      iarg++;
      if (iarg < narg) {
	dump_delay = atof(arg[iarg]);
	iarg++;
      } else error->all("Illegal dump_style command");

    } else if (strcmp(arg[iarg],"mask") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"yes") == 0) mask_flag = 1;
	else if (strcmp(arg[iarg],"no") == 0) mask_flag = 0;
	iarg++;
      } else error->all("Illegal dump_style command");
    } else break;
  }
  }

  // line = one string of concatenated keywords
  // size_one = # of keywords

  char *line = new char[MAXLINE];
  if (iarg == narg) {
    size_one = 5;
    strcpy(line,DEFAULT);
  } else {
    size_one = narg - iarg;
    line[0] = '\0';
    for (int jarg = iarg; jarg < narg; jarg++) {
      strcat(line,arg[jarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';
  }

  // allocate per-keyword memory

  vtype = new int[size_one];
  vindex = new int[size_one];
  vformat = new char*[size_one];
  pack_choice = new FnPtrPack[size_one];

  int i = 0;
  char *word = strtok(line," \0");

  while (word) {
    if (strcmp(word,"id") == 0) {
      pack_choice[i] = &Output::pack_id;
      vtype[i] = INT;
    } else if (strcmp(word,"lattice") == 0) {
      pack_choice[i] = &Output::pack_lattice;
      vtype[i] = INT;
      if (app->appclass == LATTICE && applattice->lattice == NULL)
	error->all("Dumping lattice but application does not support it");
    } else if (strcmp(word,"x") == 0) {
      pack_choice[i] = &Output::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"y") == 0) {
      pack_choice[i] = &Output::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"z") == 0) {
      pack_choice[i] = &Output::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"energy") == 0) {
      pack_choice[i] = &Output::pack_energy;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"propensity") == 0) {
      pack_choice[i] = &Output::pack_propensity;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (word[0] == 'i') {
      pack_choice[i] = &Output::pack_integer;
      vtype[i] = INT;
      vindex[i] = atoi(&word[1]);
      if (app->appclass != LATTICE)
	error->all("Invalid keyword in dump command");
      if (vindex[i] < 1 || vindex[i] > applattice->ninteger)
	error->all("Invalid keyword in dump command");
      vindex[i]--;
    } else if (word[0] == 'd') {
      pack_choice[i] = &Output::pack_double;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&word[1]) - 1;
      if (app->appclass != LATTICE)
	error->all("Invalid keyword in dump command");
      if (vindex[i] < 1 || vindex[i] > applattice->ndouble)
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

  // setup dump params

  nglobal = applattice->nglobal;
  nlocal = applattice->nlocal;
  boxxlo = applattice->boxxlo;
  boxxhi = applattice->boxxhi;
  boxylo = applattice->boxylo;
  boxyhi = applattice->boxyhi;
  boxzlo = applattice->boxzlo;
  boxzhi = applattice->boxzhi;

  if (mask_flag)
    mask = (int *) memory->smalloc(nlocal*sizeof(int),"output:mask");
}

/* ---------------------------------------------------------------------- */

void Output::add_diag(Diag *diag)
{
  ndiags++;
  diaglist = (Diag **) memory->srealloc(diaglist,ndiags*sizeof(Diag *),
					"output:diaglist");
  diaglist[ndiags-1] = diag;
}

/* ----------------------------------------------------------------------
   called only when some output is needed or when app is done
   set next output time for any output performed
   return tnext = next time any output is needed
------------------------------------------------------------------------- */

double Output::compute(double time, int done)
{
  // dump output

  if (dump_delta > 0.0 && time >= dump_time) {
    dump(time);
    dump_time = next_time(time,dump_logfreq,dump_delta,
			  dump_nrepeat,dump_scale,dump_delta);
  }

  // sflag = 1 if stats output needed

  int sflag = 0;
  if (time >= stats_time) sflag = 1;
  
  // diagnostic output, which may be driven by stats output
  
  double diag_time = app->stoptime;
  for (int i = 0; i < ndiags; i++) {
    if (diaglist[i]->stats_flag) {
      if (sflag) diaglist[i]->compute();
    } else if (time >= diaglist[i]->diag_time) {
      diaglist[i]->compute();
      diaglist[i]->diag_time = 
	next_time(time,diaglist[i]->diag_logfreq,diaglist[i]->diag_delta,
		  diaglist[i]->diag_nrepeat,diaglist[i]->diag_scale,
		  diaglist[i]->diag_delta);
      diag_time = MIN(diag_time,diaglist[i]->diag_time);
    } else diag_time = MIN(diag_time,diaglist[i]->diag_time);
  }
  
  // stats output, after diagnostics compute any needed quantities
  
  if (sflag || done) {
    stats(1);
    if (stats_delta)
      stats_time = next_time(time,stats_logfreq,stats_delta,
			     stats_nrepeat,stats_scale,stats_delta);
    else stats_time = app->stoptime;
  }

  // find next output time

  double tnext = app->stoptime;
  tnext = MIN(tnext,dump_time);
  tnext = MIN(tnext,diag_time);
  tnext = MIN(tnext,stats_time);
  return tnext;
}

/* ----------------------------------------------------------------------
   calculate next time that output of a particular kind should be performed
   return tnew = next time at which output should be done
------------------------------------------------------------------------- */

double Output::next_time(double tcurrent, int logfreq, double delta, 
			 int nrepeat, double scale, double delay)
{
  double tnew;

  if (logfreq == 0) {
    tnew = ceil(tcurrent/delta) * delta;
    if (tnew == tcurrent) tnew = tcurrent + delta;
  } else {
    double start = delta;
    while (tcurrent >= start*scale) start *= scale;
    tnew = ceil(tcurrent/start) * start;
    if (tnew == tcurrent) tnew = tcurrent + start;
    if (static_cast<int> (tnew/start) > nrepeat) tnew = start*scale;
  }

  tnew = MAX(tnew,delay);
  return tnew;
}

/* ----------------------------------------------------------------------
   print stats, including contributions from app and diagnostics
------------------------------------------------------------------------- */

void Output::stats(int timeflag)
{
  char str[2048] = {'\0'};
  char *strpnt = str;

  app->stats(strpnt);
  strpnt += strlen(strpnt);

  if (timeflag) {
    sprintf(strpnt,"%10.3g ",timer->elapsed(TIME_LOOP));
    strpnt += strlen(strpnt);
  } else {
    sprintf(strpnt,"%10.3g ",0.0);
    strpnt += strlen(strpnt);
  }

  for (int i = 0; i < ndiags; i++)
    if (diaglist[i]->stats_flag) {
      diaglist[i]->stats(strpnt);
      strpnt += strlen(strpnt);
    }

  if (me == 0) {
    if (screen)
      fprintf(screen,"%s\n",str);
    if (logfile) {
      fprintf(logfile,"%s\n",str);
      fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   print stats header, including contributions from app and diagnostics
------------------------------------------------------------------------- */

void Output::stats_header()
{
  char str[2048] = {'\0'};
  char *strpnt = str;

  app->stats_header(strpnt);
  strpnt += strlen(strpnt);

  sprintf(strpnt,"%10s","CPU");
  strpnt += strlen(strpnt);

  for (int i = 0; i < ndiags; i++)
    if (diaglist[i]->stats_flag) {
      diaglist[i]->stats_header(strpnt);
      strpnt += strlen(strpnt);
    }

  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",str);
    if (logfile) {
      fprintf(logfile,"%s\n",str);
      fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values as atom coords
------------------------------------------------------------------------- */

void Output::dump(double time)
{
  int nglobaldump,nlocaldump;

  // count sites to output with or without masking

  if (mask_flag == 0) nlocaldump = nlocal;
  else {
    for (int i = 0; i < nlocal; i++) mask[i] = 0;
    maskzeroenergy();
    nlocaldump = 0;
    for (int i = 0; i < nlocal; i++)
      if (!mask[i]) nlocaldump++;
  }

  MPI_Allreduce(&nlocaldump,&nglobaldump,1,MPI_INT,MPI_SUM,world);

  // allocate buffer for getting site info from other procs

  int nbuf = nlocaldump*size_one;
  MPI_Allreduce(&nbuf,&maxbuf,1,MPI_INT,MPI_MAX,world);
  buf = (double *) memory->smalloc(maxbuf*sizeof(double),"output:buf");

  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d %10g\n",idump,time);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nglobaldump);
    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxxlo,boxxhi);
    fprintf(fp,"%g %g\n",boxylo,boxyhi);
    fprintf(fp,"%g %g\n",boxzlo,boxzhi);
    fprintf(fp,"ITEM: ATOMS\n");
  }

  idump++;

  // pack my info into buffer

  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
  int me_size = nlocaldump*size_one;

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

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
      } else nlines = me_size/size_one;
      
      write_data(nlines,buf);
    }
    fflush(fp);

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,me_size,MPI_DOUBLE,0,0,world);
  }

  memory->sfree(buf);
}

/* ---------------------------------------------------------------------- */

void Output::write_data(int n, double *buf)
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
   one method for every keyword dump can output
   the site quantity is packed into buf starting at n with stride size_one
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Output::pack_id(int n)
{
  int *id;
  id = applattice->id;

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = id[i];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = id[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_lattice(int n)
{
  int i,j,k;
  int *lattice = applattice->lattice;

  if (mask_flag == 0) {
    for (i = 0; i < nlocal; i++) {
      buf[n] = lattice[i];
      n += size_one;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = lattice[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_x(int n)
{
  double **xyz = applattice->xyz;

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = xyz[i][0];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = xyz[i][0];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_y(int n)
{
  double **xyz = applattice->xyz;

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = xyz[i][1];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = xyz[i][1];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_z(int n)
{
  double **xyz = applattice->xyz;

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = xyz[i][2];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = xyz[i][2];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_energy(int n)
{
  int i,j,k;

  if (mask_flag == 0) {
    for (i = 0; i < nlocal; i++) {
      buf[n] = applattice->site_energy(i);
      n += size_one;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = applattice->site_energy(i);
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_propensity(int n)
{
  int i,j,k;

  if (mask_flag == 0) {
    for (i = 0; i < nlocal; i++) {
      buf[n] = applattice->site_propensity(i);
      n += size_one;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = applattice->site_propensity(i);
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_integer(int n)
{
  int *ivec = applattice->iarray[vindex[n]];

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = ivec[i];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = ivec[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_double(int n)
{
  double *dvec = applattice->darray[n];

  if (mask_flag == 0) {
    for (int i = 0; i < nlocal; i++) {
      buf[n] = dvec[i];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i]) continue;
      buf[n] = dvec[i];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Output::maskzeroenergy()
{
  for (int i = 0; i < nlocal; i++)
    mask[i] = applattice->site_energy(i) <= 0.0;
}
