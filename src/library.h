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

// prototypes for calling SPPARKS as a library

#include "mpi.h"

void spparks_open(int, char **, MPI_Comm);   // start SPPARKS w/ cmdline args
void spparks_close();                        // shut-down SPPARKS
void spparks_file(char *);                   // execute an input script
char *spparks_command(char *);               // execute a SPPARKS command
