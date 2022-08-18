/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
 Class AppPhaseFieldThinFilm - added by Fadi Abdeljawad, fabdelj@sandia.gov
 Most of this was copied from AppPotts, AppPottsNeighOnly and PottsPhaseField.
 
   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(PF/ThinFilm,AppPhaseFieldThinFilm)

#else

#ifndef SPK_APP_PHASEFIELD_THINFILM_H
#define SPK_APP_PHASEFIELD_THINFILM_H

#include "app_potts_neighonly.h"

namespace SPPARKS_NS {

class AppPhaseFieldThinFilm : public AppPottsNeighOnly {
 public:
  AppPhaseFieldThinFilm(class SPPARKS *, int, char **);
  ~AppPhaseFieldThinFilm();
  void input_app(char *, int, char **);
  void init_app();
  void grow_app();
  void site_event_rejection(int,RandomPark *);
  void app_update(double);
  double site_energy(int);
  
 protected:
  double *c;      //concentration
  double *chmPot;      //chemical potential of c
  double *mobility;      //mobility of c
  double *GB;
  double **chm_phi;
    
  int nlocal_app; //size of cnew
  int dimension;  //simulation dimension
  int *cmap;      //connectivity map for finite difference
    
  bool cmap_ready;          // flag whether the connectivity map is ready
  bool print_cmap;          // flag to print the connectivity map
  bool initialize_values;   // flag to initialize equilibrium values
  
  // constants for the free energy functional

  double  Wphi, mH, kappa, epsilon, Lphi;
  double M_o, M_vol, M_eva, M_sur, M_gb;
  int num_OP;
    
  // variables for evolving the phase field
  
  double pf_dt;           //time step for phasefield solver
  
  // methods unique to this class

  void init_values();
  void site_event_finitedifference1(int);
  void site_event_finitedifference2(int);
  double site_grain_membership(int,int);
  void setup_connectivity_map();
  void print_connectivity_map();
};

}

#endif
#endif
