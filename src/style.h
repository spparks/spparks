/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifdef AppInclude
#include "app_chemistry.h"
#include "app_ising_2d_4n.h"
#include "app_ising_2d_8n.h"
//#include "app_membrane.h"
#include "app_potts_2d_4n.h"
#include "app_potts_2d_8n.h"
#include "app_ising_3d_6n.h"
#include "app_ising_3d_26n.h"
#include "app_potts_3d_6n.h"
#include "app_potts_3d_26n.h"
#include "app_surf.h"
#include "app_test.h"
#endif

#ifdef AppClass
AppStyle(chemistry,AppChemistry)
AppStyle(ising/2d/4n,AppIsing2d4n)
//AppStyle(ising/2d/8n,AppIsing2d8n)
  //AppStyle(membrane,AppMembrane)
AppStyle(potts/2d/4n,AppPotts2d4n)
AppStyle(potts/2d/8n,AppPotts2d8n)
AppStyle(ising/3d/6n,AppIsing3d6n)
AppStyle(ising/3d/26n,AppIsing3d26n)
AppStyle(potts/3d/6n,AppPotts3d6n)
AppStyle(potts/3d/26n,AppPotts3d26n)
AppStyle(surf,AppSurf)
AppStyle(test,AppTest)
#endif

#ifdef CommandInclude
#include "shell.h"
#endif

#ifdef CommandClass
CommandStyle(shell,Shell)
#endif

#ifdef SolveInclude
#include "solve_gillespie.h"
#include "solve_next_event_alias_search.h"
#include "solve_next_event_group_search.h"
#include "solve_next_event_linear_search.h"
#include "solve_next_event_tree_search.h"
#endif

#ifdef SolveClass
SolveStyle(gillespie,SolveGillespie)
SolveStyle(next_event_alias_search,SolveNextEventAliasSearch)
SolveStyle(next_event_group_search,SolveNextEventGroupSearch)
SolveStyle(next_event_linear_search,SolveNextEventLinearSearch)
SolveStyle(next_event_tree_search,SolveNextEventTreeSearch)
#endif

#ifdef SweepInclude
#include "sweep_lattice2d.h"
#include "sweep_lattice3d.h"
#endif

#ifdef SweepClass
SweepStyle(lattice2d,SweepLattice2d)
SweepStyle(lattice3d,SweepLattice3d)
#endif

// packages

#include "style_gppt.h"

// user add-ons

#include "style_user.h"
