/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "spparks.h"
#include "input.h"

using namespace SPPARKS;

/* ----------------------------------------------------------------------
   main program to drive SPPARKS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  SPK *spk = new SPK(argc,argv,MPI_COMM_WORLD);
  spk->input->file();
  delete spk;

  MPI_Finalize();
}
