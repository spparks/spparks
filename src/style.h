
/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifdef AppInclude
#include "app_chemistry.h"
#include "app_grain.h"
#include "app_grain_strict.h"
#include "app_grain_3d.h"
#include "app_grain_3d_strict.h"
#include "app_test.h"
#endif

#ifdef AppClass
AppStyle(chemistry,AppChemistry)
AppStyle(grain,AppGrain)
AppStyle(grain_strict,AppGrainStrict)
AppStyle(grain_3d,AppGrain3D)
AppStyle(grain_3d_strict,AppGrain3DStrict)
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
#include "solve_next_event_linear_search.h"
#include "solve_next_event_tree_search.h"
#endif

#ifdef SolveClass
SolveStyle(gillespie,SolveGillespie)
SolveStyle(next_event_linear_search,SolveNextEventLinearSearch)
SolveStyle(next_event_tree_search,SolveNextEventTreeSearch)
#endif

// user add-ons

#include "style_user.h"
