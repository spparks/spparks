/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
   Class ReadTable - added by Eric Homer, ehomer@sandia.gov
   Sep 22, 2010 - This class was largely copied from the ReadSites class.

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
#include "read_table.h"
#include "table.h"
#include "app.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

#define MAXLINE 2048
#define CHUNK 1024

#define NSECTIONS 1       // change when add to header::section_keywords

/* ---------------------------------------------------------------------- */

ReadTable::ReadTable(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
}

/* ---------------------------------------------------------------------- */

ReadTable::~ReadTable()
{
  delete [] line;
  delete [] keyword;
  delete [] buffer;
}

/* ---------------------------------------------------------------------- */

void ReadTable::command(int narg, char **arg)
{
  if (app == NULL) error->all(FLERR,"Read_table command before app_style set");
  if (narg != 2) error->all(FLERR,"Illegal read_table command");
  
  char *tableName;
  
  table = (Table *) app->extract(arg[1]);
  
  if (table==NULL) error->all(FLERR,"Illegal read_table command");

  //open the file
  
  if (me == 0) {
    if (screen) fprintf(screen,"Reading table file ...\n");
    open(arg[0]);
  }

  header();
  
  // read rest of file in free format
  // if add a section keyword, add to header::section_keywords and NSECTIONS
  
  int valueflag = 0;
  
  while (strlen(keyword)) {
    if (strcmp(keyword,"Values") == 0) {
      values();
      valueflag = 1;
    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(FLERR,str);
    }
    
    parse_keyword(0);
  }
  
  // error checks
  
  if (valueflag == 0)
    error->all(FLERR,"Table file has no Values");
    
  // close file
  
  if (me == 0) {
    if (compressed) pclose(fp);
    else fclose(fp);
  }
  
  //set the flag to tell the table it's ready
  
  table->setTableReady();
}

/* ---------------------------------------------------------------------- */

void ReadTable::header()
{
  const char *section_keywords[NSECTIONS] = {"Values",};
  
  // skip 1st line of file
  
  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of table file");
  }
  
  int n,nRows=0,nCols=0,rowcolvals=0;
  char *ptr;
  
  // skip 1st line of file
  
  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of table file");
  }
  
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
    
    if (strstr(line,"size")) {
      sscanf(line,"%d %d",&nRows,&nCols);
    } 
    else if (strstr(line,"rowcolvals")) {
      rowcolvals=1;
    } else break;
  }
  
  if (nRows==0 || nCols==0)
    error->all(FLERR,"Invalid file header");
  
  // init the table now
  
  table->initTable(nRows,nCols,rowcolvals);
    
  // check that exiting string is a valid section keyword
  
  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;
  if (n == NSECTIONS) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read table values
------------------------------------------------------------------------- */

void ReadTable::values()
{
  int nRows,nCols;
  table->getTableDims(&nRows,&nCols);
  double **tableArray = table->getTablePtr();
  double *rowArray=NULL, *colArray=NULL;
  bool hasRowColVals = table->hasRowColVals();
  
  // add one to nRows and nCols if the first row and column in the file 
  // have the row column data
  
  if (hasRowColVals) {
    nCols++;
    nRows++;
    rowArray = table->getRowPtr();
    colArray = table->getColPtr();
  }
  
  int i,m,nchunk,idone;
  char *next,*buf;
  
  // read and broadcast one CHUNK of lines at a time
  // save the same table for each process
  
  int nread = 0;
  
  while (nread < nRows) {
    if (nRows-nread > CHUNK) nchunk = CHUNK;
    else nchunk = static_cast<int> (nRows - nread);
    if (me == 0) {
      char *eof;
      m = 0;
      for (i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[m],MAXLINE,fp);
        if (eof == NULL) error->one(FLERR,"Unexpected end of table file");
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
    
    if (nwords != nCols) error->all(FLERR,"Incorrect value format in table file");
    
    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      if (hasRowColVals) {
        if (nread+i  == 0) {
          double val = atof(strtok(buf," \t\n\r\f"));
          if (val != 0)
            error->all(FLERR,"First value must be zero when table has row "
                       "and column data");
          
          for (m = 0; m < nCols-1; m++) 
            colArray[m] = atof(strtok(NULL," \t\n\r\f"));
        } else {
          rowArray[nread + i-1] = atof(strtok(buf," \t\n\r\f"));
          for (m = 0; m < nCols-1; m++) 
            tableArray[nread + i-1][ m] = atof(strtok(NULL," \t\n\r\f"));
        }
      }
      else {
        tableArray[nread + i][ 0] = atof(strtok(buf," \t\n\r\f"));
          for (m = 1; m < nCols; m++) 
            tableArray[nread + i][ m] = atof(strtok(NULL," \t\n\r\f"));
      }
      
      buf = next + 1;
    }
    
    nread += nchunk;
  }
  
  if (me == 0) {
    if (screen) fprintf(screen,"  %d values\n",nRows*nCols);
    if (logfile) fprintf(logfile,"  %d values\n",nRows*nCols);
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens table file
   test if gzipped
------------------------------------------------------------------------- */

void ReadTable::open(char *file)
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
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
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

void ReadTable::parse_keyword(int first)
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
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int ReadTable::count_words(char *line)
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
