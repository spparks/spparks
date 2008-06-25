/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

// APP as a library
// C or Fortran style interface
// new APP-specific functions can be added

#include "mpi.h"
#include "library.h"
#include "app.h"
#include "input.h"

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   create an instance of SPPARKS and return pointer to it
------------------------------------------------------------------------- */

void spk_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  SPPARKS *spk = new SPPARKS(argc,argv,communicator);
  *ptr = (void *) spk;
}

/* ----------------------------------------------------------------------
   destruct an instance of SPPARKS
------------------------------------------------------------------------- */

void spk_close(void *ptr)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  delete spk;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void spk_file(void *ptr, char *str)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  spk->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *spk_command(void *ptr, char *str)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  return spk->input->one(str);
}

/* ----------------------------------------------------------------------
   add SPPARKS-specific library functions
   all must receive SPPARKS pointer as argument
------------------------------------------------------------------------- */
