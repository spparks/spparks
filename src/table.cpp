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

#include "table.h"
#include "error.h"
#include "string.h"
#include "memory.h"
#include "ctype.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Table::Table(SPPARKS *spk,char *table) : Pointers(spk) 
{
  tableName=new char[((int) strlen(table))+1];
  strcpy(tableName,table);
  
  // convert to lowercase
  
  int i=0;
  while (tableName[i]) {
    if (isupper(tableName[i]))
      tableName[i]=tolower(tableName[i]);
    i++;
  }
  
  tableReady=false;
  rowcolvals=false;
}

/* ---------------------------------------------------------------------- */

Table::~Table() 
{
  memory->destroy(table);
  delete [] tableName;
  
  if (rowcolvals) {
    memory->destroy(rowvals);
    memory->destroy(colvals);
  }
}

/* ---------------------------------------------------------------------- */

void Table::initTable(int rows,int cols,int rowcolvalsflag)
{
  nRows=rows;
  nCols=cols;
  
  memory->create(table,nRows,nCols,tableName);
  
  if (rowcolvalsflag) {
    rowcolvals = true;
    
    char *tempstr=new char[((int) strlen(tableName))+6];
    strcpy(tempstr,tableName);
    strcat(tempstr,":rows");
    memory->create(rowvals,nRows,tempstr);
    strcpy(tempstr,tableName);
    strcat(tempstr,":cols");
    memory->create(colvals,nCols,tempstr);
    
    delete [] tempstr;
  }
}

/* ---------------------------------------------------------------------- */

void Table::getTableDims(int *rows,int *cols)
{
  *rows=nRows;
  *cols=nCols;
}

/* ---------------------------------------------------------------------- */

double Table::getValue(int row,int col)
{
  if (!tableReady) error->all(FLERR,"Table values accessed before table initialized");
  return table[row][col];
}

/* ---------------------------------------------------------------------- */

bool Table::tableNameCompare(const char *table)
{
  char * tableLower=new char[((int) strlen(table))+1];
  strcpy(tableLower,table);
  
  // convert to lowercase
  
  int i=0;
  while (tableLower[i]) {
    if (isupper(tableLower[i]))
      tableLower[i]=tolower(tableLower[i]);
    i++;
  }
  
  bool cmp;
  if (strcmp(tableLower,tableName) == 0) cmp=true;
  else cmp=false;
  
  delete tableLower;
  return cmp;
}

/* ---------------------------------------------------------------------- */

void Table::printTable()
{
  if (rowcolvals) {
    printf("         ");
    for (int j=0; j<nCols; j++)
      printf("%.5f  ",colvals[j]);
    printf("\n");
  }
  for (int i=0; i<nRows; i++) {
    if (rowcolvals) printf("%.5f  ",rowvals[i]);
    for (int j=0; j<nCols; j++) {
      printf("%.5f  ",table[i][j]);
    }
    printf("\n");
  }
}

/* ---------------------------------------------------------------------- */

void Table::setTableReady()
{
  tableReady=true;
}
