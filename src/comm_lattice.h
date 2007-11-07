/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef COMM_LATTICE_H
#define COMM_LATTICE_H

#include "sysptr.h"

namespace SPPARKS {

class CommLattice : protected SysPtr {
 public:
  CommLattice(class SPK *);
  ~CommLattice();
  void init(const int,
	    const int, const int, const int, const int,
	    const int, const int);
  void init(const int,
	    const int, const int, const int, const int, const int, const int,
	    const int, const int);
  void all(int *);
  void sector(int *, const int);
  void reverse_sector(int *, const int);

 private:
  int me,nprocs;
  int dimension;
  int nlocal;
  int procwest,proceast,procsouth,procnorth,procdown,procup;
  int delghost,dellocal; // thickness of ghost and local communication layers 

  void setup_swapinfo();
  void setup_reverseinfo();
  void sector_multilayer(int *, const int);
  void sector_multilayer_destroy(int *, const int);
  void reverse_sector_multilayer(int *, const int);
  void all_multilayer(int *);
  void allocate_swap(const int, const int);          // allocate swap arrays
  void free_swap();                                  // free swap arrays
};

}

#endif

