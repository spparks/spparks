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

#ifdef APP_CLASS

AppStyle(potts/pfm,AppPottsPhaseField)

#else

#ifndef SPK_APP_POTTS_PHASEFIELD_H
#define SPK_APP_POTTS_PHASEFIELD_H

#include "app_potts_neighonly.h"

namespace SPPARKS_NS {

class AppPottsPhaseField : public AppPottsNeighOnly {
 public:
  
  AppPottsPhaseField(class SPPARKS *, int, char **);
  ~AppPottsPhaseField();

  void init_app();
  void setup_end_app();
  void grow_app();

  void site_event_rejection(int,RandomPark *);
  void app_update(double);
  double site_energy(int);
  
 protected:
  
  double *c;      //concentration
  double *cnew;   //concentration following evolution - kept locally
  int nlocal_app; //size of cnew
  int dimension;  //simulation dimension
  int *cmap;      //connectivity map for finite difference
  
  int phaseChangeInt; //integer equal to half the # of spins to distinguish phases
  int *phase;         //pointer to phase variable
  
  bool cmap_ready;                // flag to indicate whether the connectivity map is ready
  bool print_cmap;                // flag to print the connectivity map
  bool enforceConcentrationLimits;// flag to force concentration to stay between 0 and 1
  bool pf_resetfield;             // flag to reset concentration (constant sink and source)
  bool initialize_values;         // flag to indicate whether to initialize equilibrium values
  
  //constants for the free energy functional
  double gamma,kappaC,M_c;
  double c_1,c_2,c_3,c_4;
  double a_1,a_2;
  
  //variables for evolving the phase field
  double dt_phasefield;           //time step for phasefield solver
  double dt_phasefield_mult;      // number of PF iterations per rkmc event
  
  //variables for reseting the concentration (constant sink and source)
  int *pf_resetlist;
  double *pf_resetlistvals;
  int pf_nresetlist;
  
  //variables to warn when concentration has deviated outside [0,1]
  int warn_concentration_deviation;     //local proc flag
  int warn_concentration_deviation_all; //global flag
  
  
  // methods unique to this class
  void init_values();
  void site_event_finitedifference(int);
  double site_energy_no_gradient(int i);
  void check_phasefield_stability(int,double);
  
  double site_grain_membership(int,int);
  void set_site_phase(int);
  
  void setup_connectivity_map();
  void print_connectivity_map();
  
  void set_phasefield_resetfield();
};

}

#endif
#endif
