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

/* 
   C or Fortran style library interface to SPPARKS
   new SPPARKS-specific functions can be added
*/

#include "mpi.h"

void spparks_open(int, char **, MPI_Comm, void **);  /* start-up SPPARKS */
void spparks_close(void *);                          /* shut-down SPPARKS */
void spparks_file(void *, char *);                   /* run an input script */
char *spparks_command(void *, char *);               /* execute a command */

void *spparks_extract(void *, char *);
double spparks_energy(void *);
