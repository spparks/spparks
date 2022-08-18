/* ----------------------------------------------------------------------
 SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
 http://www.cs.sandia.gov/~sjplimp/spparks.html
 Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
 Eric Homer
 
 Copyright (2008) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under 
 the GNU General Public License.
 
 See the README file in the top-level SPPARKS directory.
 ------------------------------------------------------------------------- */

#ifndef SPK_TABLE_H
#define SPK_TABLE_H

#include "pointers.h"
#include "string.h"


namespace SPPARKS_NS {
  
class Table : protected Pointers {
  public:
    Table(class SPPARKS *,char *);
    ~Table();
    void initTable(int,int,int);
    double** getTablePtr() {return table;}
    double* getRowPtr() {return rowvals;}
    double* getColPtr() {return colvals;}
    void getTableDims(int*,int*);
    double getValue(int,int);
    bool tableNameCompare(const char *);
    void printTable();
    void setTableReady();
    bool getTableReady() {return tableReady;}
    bool hasRowColVals() {return rowcolvals;}
  private:
    double** table;
    double *rowvals, *colvals;
    bool rowcolvals;
    int nRows,nCols;
    char* tableName;
    bool tableReady;
  };
  
}

#endif