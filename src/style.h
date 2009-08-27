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

#ifdef AppInclude
#include "app_chemistry.h"
#include "app_diffusion.h"
#include "app_ising.h"
#include "app_ising_single.h"
#include "app_membrane.h"
#include "app_potts.h"
#include "app_potts_neigh.h"
#include "app_potts_neighonly.h"
#include "app_potts_pin.h"
#include "app_potts_variable.h"
#include "app_test_group.h"
#endif

#ifdef AppClass
AppStyle(chemistry,AppChemistry)
AppStyle(diffusion,AppDiffusion)
AppStyle(ising,AppIsing)
AppStyle(ising/single,AppIsingSingle)
AppStyle(membrane,AppMembrane)
AppStyle(potts,AppPotts)
AppStyle(potts/neigh,AppPottsNeigh)
AppStyle(potts/neighonly,AppPottsNeighOnly)
AppStyle(potts/pin,AppPottsPin)
AppStyle(potts/variable,AppPottsVariable)
AppStyle(test/group,AppTestGroup)
#endif

#ifdef CommandInclude
#include "shell.h"
#endif

#ifdef CommandClass
CommandStyle(shell,Shell)
#endif

#ifdef SolveInclude
#include "solve_group.h"
#include "solve_linear.h"
#include "solve_tree.h"
#endif

#ifdef SolveClass
SolveStyle(group,SolveGroup)
SolveStyle(linear,SolveLinear)
SolveStyle(tree,SolveTree)
#endif

#ifdef DiagInclude
#include "diag_cluster.h"
#include "diag_diffusion.h"
#include "diag_energy.h"
#include "diag_propensity.h"
#endif

#ifdef DiagClass
DiagStyle(cluster,DiagCluster)
DiagStyle(diffusion,DiagDiffusion)
DiagStyle(energy,DiagEnergy)
DiagStyle(propensity,DiagPropensity)
#endif

// packages

//#include "style_gppt.h"

// user add-ons

#include "style_user.h"
