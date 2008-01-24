/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifdef AppInclude
#include "app_chemistry.h"
#include "app_ising.h"
#include "app_ising_2d_4n.h"
#include "app_ising_2d_8n.h"
#include "app_ising_3d_6n.h"
#include "app_ising_3d_26n.h"
#include "app_membrane.h"
#include "app_membrane_2d.h"
#include "app_potts.h"
#include "app_potts_variable.h"
#include "app_potts_2d_4n.h"
#include "app_potts_2d_8n.h"
#include "app_potts_2d_24n.h"
#include "app_potts_3d_6n.h"
#include "app_potts_3d_12n.h"
#include "app_potts_3d_26n.h"
#include "app_surf.h"
#include "app_template.h"
#include "app_test_group.h"
#include "app_test2_group.h"
#endif

#ifdef AppClass
AppStyle(chemistry,AppChemistry)
AppStyle(ising,AppIsing)
AppStyle(ising/2d/4n,AppIsing2d4n)
AppStyle(ising/2d/8n,AppIsing2d8n)
AppStyle(ising/3d/6n,AppIsing3d6n)
AppStyle(ising/3d/26n,AppIsing3d26n)
AppStyle(membrane,AppMembrane)
AppStyle(membrane/2d,AppMembrane2d)
AppStyle(potts,AppPotts)
AppStyle(potts/variable,AppPottsVariable)
AppStyle(potts/2d/4n,AppPotts2d4n)
AppStyle(potts/2d/8n,AppPotts2d8n)
AppStyle(potts/2d/24n,AppPotts2d24n)
AppStyle(potts/3d/6n,AppPotts3d6n)
AppStyle(potts/3d/12n,AppPotts3d12n)
AppStyle(potts/3d/26n,AppPotts3d26n)
AppStyle(surf,AppSurf)
AppStyle(template,AppTemplate)
AppStyle(test/group,AppTestGroup)
AppStyle(test2/group,AppTest2Group)
#endif

#ifdef CommandInclude
#include "shell.h"
#endif

#ifdef CommandClass
CommandStyle(shell,Shell)
#endif

#ifdef SolveInclude
#include "solve_alias.h"
#include "solve_group.h"
#include "solve_group2.h"
#include "solve_linear.h"
#include "solve_tree.h"
#endif

#ifdef SolveClass
SolveStyle(alias,SolveAlias)
SolveStyle(group,SolveGroup)
SolveStyle(group2,SolveGroup2)
SolveStyle(linear,SolveLinear)
SolveStyle(tree,SolveTree)
#endif

#ifdef SweepInclude
#include "sweep_lattice.h"
#include "sweep_lattice2d.h"
#include "sweep_lattice3d.h"
#endif

#ifdef SweepClass
SweepStyle(lattice,SweepLattice)
SweepStyle(lattice2d,SweepLattice2d)
SweepStyle(lattice3d,SweepLattice3d)
#endif

#ifdef DiagInclude
#include "diag_cluster.h"
#include "diag_cluster2d.h"
#include "diag_cluster3d.h"
#include "diag_eprof3d.h"
#include "diag_energy.h"
#include "diag_energy2d.h"
#include "diag_energy3d.h"
#endif

#ifdef DiagClass
DiagStyle(cluster,DiagCluster)
DiagStyle(cluster2d,DiagCluster2d)
DiagStyle(cluster3d,DiagCluster3d)
DiagStyle(energy,DiagEnergy)
DiagStyle(energy2d,DiagEnergy2d)
DiagStyle(energy3d,DiagEnergy3d)
DiagStyle(eprof3d,DiagEprof3d)
#endif

// packages

#include "style_gppt.h"

// user add-ons

#include "style_user.h"
