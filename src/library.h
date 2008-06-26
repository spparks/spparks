/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

// prototypes for calling SPPARKS as a library

#include "mpi.h"

void spparks_open(int, char **, MPI_Comm);   // start SPPARKS w/ cmdline args
void spparks_close();                        // shut-down SPPARKS
void spparks_file(char *);                   // execute an input script
char *spparks_command(char *);               // execute a SPPARKS command
