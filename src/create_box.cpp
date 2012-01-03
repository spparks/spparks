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

#include "stdlib.h"
#include "string.h"
#include "create_box.h"
#include "app.h"
#include "domain.h"
#include "region.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(SPPARKS *spk) : Pointers(spk) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (app == NULL) error->all(FLERR,"Create_box command before app_style set");

  if (narg != 1) error->all(FLERR,"Illegal create_box command");

  if (app->appclass == App::GENERAL)
    error->all(FLERR,"Cannot create box with this application style");
  if (domain->box_exist) 
    error->all(FLERR,"Cannot create box after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (domain->dimension == 1 && 
      (domain->yperiodic == 0 || domain->zperiodic == 0))
    error->all(FLERR,"Cannot run 1d simulation with nonperiodic Y or Z dimension");

  // region check

  int iregion = domain->find_region(arg[0]);
  if (iregion == -1) error->all(FLERR,"Create_box region ID does not exist");
  if (domain->regions[iregion]->interior == 0)
    error->all(FLERR,"Create_box region must be of type inside");

  // setup simulation box from region extent

  domain->boxxlo = domain->regions[iregion]->extent_xlo;
  domain->boxxhi = domain->regions[iregion]->extent_xhi;
  domain->boxylo = domain->regions[iregion]->extent_ylo;
  domain->boxyhi = domain->regions[iregion]->extent_yhi;
  domain->boxzlo = domain->regions[iregion]->extent_zlo;
  domain->boxzhi = domain->regions[iregion]->extent_zhi;

  domain->set_box();
  domain->box_exist = 1;

  if (domain->me == 0)
    if (screen) {
      if (screen) fprintf(screen,"Created box = (%g %g %g) to (%g %g %g)\n",
			  domain->boxxlo,domain->boxylo,domain->boxzlo,
			  domain->boxxhi,domain->boxyhi,domain->boxzhi);
      if (logfile) fprintf(logfile,"Created box = (%g %g %g) to (%g %g %g)\n",
			   domain->boxxlo,domain->boxylo,domain->boxzlo,
			   domain->boxxhi,domain->boxyhi,domain->boxzhi);
    }

  if (domain->dimension == 1) domain->procs2domain_1d();
  if (domain->dimension == 2) domain->procs2domain_2d();
  if (domain->dimension == 3) domain->procs2domain_3d();
  
  if (domain->me == 0)
    if (screen) {
      if (screen) fprintf(screen,"  %d by %d by %d processor grid\n",
			  domain->procgrid[0],domain->procgrid[1],
			  domain->procgrid[2]);
      if (logfile) fprintf(logfile,"  %d by %d by %d processor grid\n",
			   domain->procgrid[0],domain->procgrid[1],
			   domain->procgrid[2]);
    }
}
