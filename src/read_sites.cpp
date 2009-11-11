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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_sites.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

#include <map>

using namespace SPPARKS_NS;

#define MAXLINE 256
#define CHUNK 1024
#define DELTA 4

#define NSECTIONS 3       // change when add to header::section_keywords

#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ReadSites::ReadSites(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;
}

/* ---------------------------------------------------------------------- */

ReadSites::~ReadSites()
{
  delete [] line;
  delete [] keyword;
  delete [] buffer;
  memory->sfree(arg);
}

/* ---------------------------------------------------------------------- */

void ReadSites::command(int narg, char **arg)
{
  if (app == NULL) error->all("Read_sites command before app_style set");

  if (narg != 1) error->all("Illegal read_sites command");

  // read header info

  if (me == 0) {
    if (screen) fprintf(screen,"Reading site file ...\n");
    open(arg[0]);
  }
  header();

  // if box already exists, test that data file box is consistent
  // else create it

  if (domain->box_exist) {
    if (fabs(domain->boxxlo-boxxlo) > EPSILON ||
	fabs(domain->boxylo-boxylo) > EPSILON ||
	fabs(domain->boxzlo-boxzlo) > EPSILON ||
	fabs(domain->boxxhi-boxxhi) > EPSILON ||
	fabs(domain->boxyhi-boxyhi) > EPSILON ||
	fabs(domain->boxzhi-boxzhi) > EPSILON)
      error->all("Read_sites simluation box different that current box");

  } else {
    domain->boxxlo = boxxlo;
    domain->boxylo = boxylo;
    domain->boxzlo = boxzlo;
    domain->boxxhi = boxxhi;
    domain->boxyhi = boxyhi;
    domain->boxzhi = boxzhi;

    domain->set_box();
    domain->box_exist = 1;
  }

  // if sites already exist, test that data file nglobal is consistent

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  }

  if (app->sites_exist) {
    if (latticeflag) {
      if (nglobal != applattice->nglobal)
	error->all("Number of sites does not match existing sites");
    } else {
      if (nglobal != appoff->nglobal)
	error->all("Number of sites does not match existing sites");
    }
  } else app->sites_exist = 1;

  // read rest of file in free format
  // if add a section keyword, add to header::section_keywords and NSECTIONS

  int sitesflag = 0;

  while (strlen(keyword)) {
    if (strcmp(keyword,"Sites") == 0) {
      if (app->sites_exist) 
	error->all("Cannot read sites after sites already exist");
      sites();
      sitesflag = 1;
      app->sites_exist = 1;

    } else if (strcmp(keyword,"Neighbors") == 0) {
      if (sitesflag == 0) error->all("Must read Sites before Neighbors");
      if (latticeflag == 0) 
	error->all("Can only read neighbors for on-lattice applications");
      if (maxneigh <= 0) 
	error->all("Cannot read neighbors unless max neighbors is set");

      applattice->maxneigh = maxneigh;
      applattice->grow(app->nlocal);
      neighbors();

    } else if (strcmp(keyword,"Values") == 0) {
      if (app->sites_exist == 0 && sitesflag == 0) 
	error->all("Must read Sites before Values");
      values();

    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(str);
    }

    parse_keyword(0);
  }

  // close file

  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   if flag = 0, only called by proc 0
   if flag = 1, called by all procs so bcast lines as read them
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void ReadSites::header()
{
  int n;
  char *ptr;

  char *section_keywords[NSECTIONS] = {"Sites","Neighbors","Values",};
  
  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one("Unexpected end of data file");
  }

  // defaults

  nglobal = 0;
  maxneigh = 0;
  boxxlo = boxylo = boxzlo = -0.5;
  boxxhi = boxyhi = boxzhi = 0.5;

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    // bcast line

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable

    if (strstr(line,"sites")) sscanf(line,"%d",&nglobal);
    else if (strstr(line,"max neighbors")) sscanf(line,"%d",&maxneigh);
    else if (strstr(line,"xlo xhi")) sscanf(line,"%lg %lg",&boxxlo,&boxxhi);
    else if (strstr(line,"ylo yhi")) sscanf(line,"%lg %lg",&boxylo,&boxyhi);
    else if (strstr(line,"zlo zhi")) sscanf(line,"%lg %lg",&boxzlo,&boxzhi);
    else break;
  }

  // check that exiting string is a valid section keyword

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;
  if (n == NSECTIONS) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->all(str);
  }
}

/* ----------------------------------------------------------------------
   read all sites
   accumulate nread in double precision to allow nglobal > 2^31
------------------------------------------------------------------------- */

void ReadSites::sites()
{
  int i,m,nchunk,id;
  double x,y,z;
  char *values[4];
  char *next,*buf;

  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;

  // read and broadcast one CHUNK of lines at a time
  // add a site if I own its coords

  int nread = 0;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = static_cast<int> (nglobal - nread);
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of data file");
	m += strlen(&buffer[m]);
      }
      buffer[m++] = '\n';
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);
    buf = buffer;

    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = count_words(buf);
    *next = '\n';

    if (nwords != 4) error->all("Incorrect site format in data file");

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');

      values[0] = strtok(buf," \t\n\r\f");
      values[1] = strtok(NULL," \t\n\r\f");
      values[2] = strtok(NULL," \t\n\r\f");
      values[3] = strtok(NULL," \t\n\r\f");

      id = atoi(values[0]);
      x = atof(values[1]);
      y = atof(values[2]);
      z = atof(values[3]);
      
      if (x >= subxlo && x < subxhi &&
	  y >= subylo && y < subyhi &&
	  z >= subzlo && z < subzhi) {
	if (latticeflag) applattice->add_site(id,x,y,z);
	else appoff->add_site(id,x,y,z);
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  // check that all sites were assigned correctly

  MPI_Allreduce(&app->nlocal,&app->nglobal,1,MPI_INT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d sites\n",app->nglobal);
    if (logfile) fprintf(logfile,"  %d sites\n",app->nglobal);
  }

  if (nglobal != app->nglobal) 
    error->all("Did not assign all sites correctly");
  
  // check that sites IDs range from 1 to nglobal
  // not checking if site IDs are unique
  
  int flag = 0;
  for (int i = 0; i < app->nlocal; i++)
    if (app->id[i] <= 0 || app->id[i] > nglobal) flag = 1;
  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all("Invalid site ID in Sites section of data file");
}

/* ----------------------------------------------------------------------
   read all neighbors of sites
   to find atoms, must build atom map if not a molecular system 
   accumulate nread in double precision to allow natoms > 2^31
------------------------------------------------------------------------- */

void ReadSites::neighbors()
{
  int i,m,nchunk,idone,nvalues;
  char *next,*buf;

  // put my entire list of owned site IDs in a hash

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;

  int *id = app->id;
  int nlocal = app->nlocal;

  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  // read and broadcast one CHUNK of lines at a time
  // store site's neighbors if I own its ID

  char **values = new char*[maxneigh+1];

  int nread = 0;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = static_cast<int> (nglobal - nread);
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of data file");
	m += strlen(&buffer[m]);
      }
      buffer[m++] = '\n';
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);

    buf = buffer;
    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');

      idone = atoi(strtok(buf," \t\n\r\f"));
      loc = hash.find(idone);

      if (loc != hash.end()) {
	nvalues = 0;
	while (nvalues <= maxneigh) {
	  values[nvalues] = strtok(NULL," \t\n\r\f");
	  if (values[nvalues] == NULL) break;
	  nvalues++;
	}

	if (nvalues > maxneigh) error->one("Too many neighbors per site");
	applattice->add_neighbors(loc->second,nvalues,values);
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  delete [] values;

  int ncount = 0;
  for (i = 0; i < app->nlocal; i++)
    ncount += applattice->numneigh[i];
  int ntotal;
  MPI_Allreduce(&ncount,&ntotal,1,MPI_INT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d neighbors\n",ntotal);
    if (logfile) fprintf(logfile,"  %d neighbors\n",ntotal);
  }
}

/* ----------------------------------------------------------------------
   read all per-site values
   to find atoms, must build atom map if not a molecular system 
   accumulate nread in double precision to allow natoms > 2^31
------------------------------------------------------------------------- */

void ReadSites::values()
{
  int i,m,nchunk,idone;
  char *next,*buf;

  // put my entire list of owned site IDs in a hash

  std::map<int,int>::iterator loc;
  std::map<int,int> hash;

  int *id = app->id;
  int nlocal = app->nlocal;

  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));

  // read and broadcast one CHUNK of lines at a time
  // store site's values if I own its ID

  int nvalues = app->ninteger + app->ndouble;
  char **values = new char*[nvalues];

  int nread = 0;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = static_cast<int> (nglobal - nread);
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
	eof = fgets(&buffer[m],MAXLINE,fp);
	if (eof == NULL) error->one("Unexpected end of data file");
	m += strlen(&buffer[m]);
      }
      buffer[m++] = '\n';
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_CHAR,0,world);
    buf = buffer;

    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = count_words(buf);
    *next = '\n';

    if (nwords != nvalues+1) error->all("Incorrect value format in data file");

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');

      idone = atoi(strtok(buf," \t\n\r\f"));
      loc = hash.find(idone);

      if (loc != hash.end()) {
	for (m = 0; m < nvalues; m++) values[m] = strtok(NULL," \t\n\r\f");
	if (latticeflag) applattice->add_values(loc->second,values);
	else appoff->add_values(loc->second,values);
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  delete [] values;

  if (me == 0) {
    if (screen) fprintf(screen,"  %d values\n",nglobal*nvalues);
    if (logfile) fprintf(logfile,"  %d values\n",nglobal*nvalues);
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void ReadSites::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    char gunzip[128];
    sprintf(gunzip,"gunzip -c %s",file);
    fp = popen(gunzip,"r");
#else
    error->one("Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(str);
  }
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void ReadSites::parse_keyword(int first)
{
  int eof = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t' 
	 || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   parse a line of coeffs into words, storing them in narg,arg
   trim anything from '#' onward
   word strings remain in line, are not copied
   if addflag, duplicate 1st word, so pair_coeff "2" looks like "2 2"
------------------------------------------------------------------------- */

void ReadSites::parse_coeffs(int addflag, char *line)
{
  char *ptr;
  if (ptr = strchr(line,'#')) *ptr = '\0';

  narg = 0;
  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) 
	memory->srealloc(arg,maxarg*sizeof(char *),"read_sites:arg");
    }
    arg[narg++] = word;
    if (addflag && narg == 1) continue;
    word = strtok(NULL," \t\n\r\f");
  }
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int ReadSites::count_words(char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}
