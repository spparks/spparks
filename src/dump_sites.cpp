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

#include "dump_sites.h"
#include "app.h"
#include "domain.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // in other dump files

/* ---------------------------------------------------------------------- */

DumpSites::DumpSites(SPPARKS *spk, int narg, char **arg) : 
  DumpText(spk, narg, arg)
{
  // error check for allowed dump sites output

  if (!multifile)
    error->all(FLERR,"Dump sites must write one file per snapshot");
  if (multiproc) 
    error->all(FLERR,"Dump sites cannot write multiple files per snapshot");
  if (binary) error->all(FLERR,"Dump sites cannot write to binary files");

  // require ID first, then one or more per-site values

  if (fields[0] != ID) 
    error->all(FLERR,"Dump sites requires ID as first field");
  if (size_one < 2) error->all(FLERR,"Dump sites requires two or more fields");


}

/* ---------------------------------------------------------------------- */

void DumpSites::write_header(bigint ndump, double time)
{
  if (ndump != app->nglobal)
    error->all(FLERR,"Dump sites must output info for all sites");

  if (me) return;

  fprintf(fp,"Site file written by dump sites %s command at time: %d %g\n\n",
          id,idump,time);
  fprintf(fp,"%d dimension\n",domain->dimension);
  fprintf(fp,BIGINT_FORMAT " sites\n",ndump);
  fprintf(fp,"%s values\n",columns_orig);
  fprintf(fp,"%g %g xlo xhi\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g ylo yhi\n",boxylo,boxyhi);
  fprintf(fp,"%g %g zlo zhi\n",boxzlo,boxzhi);
  fprintf(fp,"\nValues\n\n");
}
