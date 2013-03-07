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

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void spparks_open(int, char **, MPI_Comm, void **);
void spparks_open_no_mpi(int, char **, void **);
void spparks_close(void *);
void spparks_file(void *, char *);
char *spparks_command(void *, char *);

void *spparks_extract(void *, char *);
double spparks_energy(void *);

#ifdef __cplusplus
}
#endif
