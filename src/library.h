/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

// prototypes for calling SPPARKS as a library

#include "mpi.h"

void spk_open(int, char **, MPI_Comm); // start SPPARKS w/ command-line args
void spk_close();                      // shut-down SPPARKS
void spk_file(char *);                 // execute an input script
char *spk_command(char *);             // execute a single SPPARKS command
