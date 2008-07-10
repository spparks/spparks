/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifdef AppInclude
#include "app_chemistry.h"
#include "app_ising.h"
#include "app_ising_exchange.h"
#include "app_ising_2d_4n.h"
#include "app_ising_2d_4n_exchange.h"
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
#include "app_test_group.h"
#endif

#ifdef AppClass
AppStyle(chemistry,AppChemistry)
AppStyle(ising,AppIsing)
AppStyle(ising/exchange,AppIsingExchange)
AppStyle(ising/2d/4n,AppIsing2d4n)
AppStyle(ising/2d/4n/exchange,AppIsing2d4nExchange)
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

//#include "style_gppt.h"

// user add-ons

#include "style_user.h"
