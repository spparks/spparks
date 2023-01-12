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
AppStyle(phasefield/potts/table,AppPhaseFieldPottsTable)

#else

#ifndef SPK_APP_PHASEFIELD_POTTS_TABLE_H
#define SPK_APP_PHASEFIELD_POTTS_TABLE_H

#include "app_potts_neighonly.h"

namespace SPPARKS_NS {
  
 class AppPhaseFieldPottsTable : public AppPottsNeighOnly {
  public:
  
  AppPhaseFieldPottsTable(class SPPARKS *, int, char **);
  ~AppPhaseFieldPottsTable();

  void input_app(char *, int, char **);
  void init_app();
  void setup_end_app();
  void grow_app();
  void *extract_app(char *);

  void site_event_rejection(int,RandomPark *);
  void app_update(double);
  double site_energy(int);
  
 protected:
  
  double *c;      //concentration
  double *cnew;   //concentration following evolution - kept locally
  double *T;      //temperature
  double *Tnew;   //temperature following evolution - kept locally
  
  int nlocal_app; //size of cnew
  int dimension;  //simulation dimension
  int *cmap;      //connectivity map for finite difference
  double latconst;  //set as the lattice constant from lattice
  
  int *phase;         //pointer to phase variable
  
  bool cmap_ready;                // flag to indicate whether the connectivity map is ready
  bool print_cmap;                // flag to print the connectivity map
  bool init_site_phase;
  bool enforceConcentrationLimits;// flag to force concentration to stay between 0 and 1
  bool pf_resetfield;             // flag to reset concentration (constant sink and source)
  bool run_nucleation;            // flag to indicate whether to run nuclation
  double scale_potts_energy;      // value by which to scale the potts interface energy
  double scale_free_energy;       // value by which to scale the free energy
  double scale_comp_energy_penalty; // value by which to scale the composition energy penalty when it is out of bounds
  double nucleation_rate;         // nucleation rate
  
  //constants for the Cahn-Hilliard equation
  double mobility,ch_energy;//mobililty and Cahn-Hilliard factor for the energy
  
  //constants for heat and mass transfer
  double thermal_conductivity;
  double heat_of_transport;
  double specific_heat;
  
  
  struct Phase {
    char *id;
    int phase;
    int spin_start,spin_end,spin_range;
    class Table *table;
    int nRow,nCol;
    double *rowvals,*colvals;
    double dC,dT;
    double **G, **dGdC;
    double scaleval;
  };
  
  int nphases;
  Phase **phases;           // list of phases
  
  //variables for evolving the phase field
  double dt_phasefield;           //time step for phasefield solver
  double dt_phasefield_mult;      // number of PF iterations per rkmc event
  
  //variables for reseting the concentration (constant sink and source)
  int *pf_resetlist;
  double *pf_resetlistvals;
  int pf_nresetlist,pf_maxresetlist;
  
  //variables to warn when concentration has deviated outside [0,1]
  int warn_concentration_deviation;     //local proc flag
  int warn_concentration_deviation_all; //global flag  
  
  //variables to warn when temperature has attemped to go below 0
  int warn_temperature_deviation;     //local proc flag
  int warn_temperature_deviation_all; //global flag
  
  // methods unique to this class
  void site_event_finitedifference(int);
  double site_energy_without_phasefield(int i);
  double site_energy_bulk(int i);
  double site_energy_potts(int i);
  double site_energy_phasefield(int i);
  
  int find_phase(char *);
  void set_site_phase(int);
  void process_phases();
  
  void setup_connectivity_map();
  void print_connectivity_map();
  
  void set_pf_resetlist(int, char **);
  double bilinear_interp(int, int);
  int binarySearch(double *, int, int, double);
};

}

#endif
#endif
