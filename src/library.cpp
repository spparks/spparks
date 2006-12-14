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

using namespace SPPARKS;

/* ----------------------------------------------------------------------
   create an instance of SPK and return pointer to it
------------------------------------------------------------------------- */

void spk_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  SPK *spk = new SPK(argc,argv,communicator);
  *ptr = (void *) spk;
}

/* ----------------------------------------------------------------------
   destruct an instance of SPK
------------------------------------------------------------------------- */

void spk_close(void *ptr)
{
  SPK *spk = (SPK *) ptr;
  delete spk;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void spk_file(void *ptr, char *str)
{
  SPK *spk = (SPK *) ptr;
  spk->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *spk_command(void *ptr, char *str)
{
  SPK *spk = (SPK *) ptr;
  return spk->input->one(str);
}

/* ----------------------------------------------------------------------
   add SPK-specific library functions
   all must receive SPK pointer as argument
------------------------------------------------------------------------- */
