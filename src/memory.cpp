/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   Website
   https://spparks.github.io/

   See authors 
   https://spparks.github.io/authors.html

   Copyright(C) 1999-2025 National Technology & Engineering Solutions
   of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.  See 
   LICENSE in top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "spktype.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(SPPARKS *spk) : Pointers(spk) {}

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate " BIGINT_FORMAT " bytes for array %s",
	    nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc 
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return NULL;
  }

  ptr = realloc(ptr,nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate " BIGINT_FORMAT " bytes for array %s",
	    nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str,"Cannot create/grow a vector/array of pointers for %s",name);
  error->one(FLERR,str);
}
