
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
AppStyle(sinter,AppSinter)

#else

#ifndef APP_SINTER_H
#define APP_SINTER_H

#include <map>
#include <vector>
using std::vector;

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppSinter : public AppLattice {
	friend class DiagSinterDensity;
	friend class DiagSinterDensityTK;
	friend class DiagSinterFreeEnergy;
	friend class DiagSinterFreeEnergyPore;
	
	enum {FRAME=-1, VACANT};
			
 public:
  AppSinter(class SPPARKS *, int, char **);
  ~AppSinter();
  
  // Application Initialization
  void grow_app();
  void input_app(char *command, int narg, char **arg);
  void init_app();
  void setup_app();
  void initialize_parameters_density_calculation();
  void overimpose_frame();
  void check_space_initialization();  
  
  // Functions required as class derives from AppLattice
  double site_energy(int);
  double site_propensity(int);					// Not used (part of rejection-free KMC)
  void site_event(int, class RandomPark *);	// Not used (part of rejection-free KMC)
  void site_event_rejection(int, class RandomPark *);

  // Change introduced in AppLattice to be able to stop and synchronize annihilations after sweepping sector
  void iterate_rejection(double stoptime);


  // Mapping between real space, simulation grid and processor grid
  // NOTE : Mapping to processor depends on the actual implementation of AppLattice and Domain classes 
  //        --> review if changes to these base clases are made
  int xyz_to_processor_and_grid( double x, double y, double z, int &grid_x, int &grid_y, int &grid_z );
  int xyz_to_processor_and_global( double x, double y, double z, int &iglobal );
  int xyz_to_local( double x, double y, double z );
  void global_to_grid( int iglobal, int & i, int & j, int & k );

  // Function to select spin from neighboring grains minimizing energy
  void choose_neighbor_grain_site_minimizing_energy( int i, int neighstate, RandomPark *random, double & min_gs_energy, int & newstate );

  // Functions for vacancies and annihilations
  void make_vacancy(int i, RandomPark *random);		// Serial / Parallel
  void annihilate(int i, RandomPark *random);		// Serial / Parallel
  void register_annihilation_event(int i);			// Parallel
  void update_annihilations(RandomPark *random);	// Parallel

  // Processing annihilations SERIAL version
  bool vacancy_adjacent_grain( int i, RandomPark *random, int & adjacent_grain_spin, int & adjacent_grain_start );	// Serial / Parallel
  double calculate_mass_center_adjacent_grain(int i, int isite, vector<double> & cm );
  void calculate_vacancy_new_position( int i, vector<double> & cm, vector<double> & p, double &mint, vector<double> & new_pos ); // Serial / Parallel
  void collapse_in_direction( int ind_vac, vector<double> & p, double limit_t, vector<double> & new_pos );
  
  // Processing annihilations PARALLEL version
  void generate_annihilation_spin_list(RandomPark *random);
  void register_collapsing_event(int ind_vac, vector<double> & p, double limit_t, vector<double> & new_pos);
//  void calculate_mass_center_distributed_grains( const int SIZE_LIST, const int SIZE_GROUP, vector<double> & cm_grains );
//  void generate_local_vacancy_collapsing_list( const int SIZE_LIST, const int SIZE_GROUP, vector<double> & cm_grains );
  void generate_local_vacancy_collapsing_list( const int SIZE_LIST );
  void filter_local_vacancy_collapsing_list();
  void collapse_pending_list ( const int SIZE_LIST ); 
 
  // Functions for grain boundary calculations
  double calculate_gb_update(double current_time);			  
  double calculate_gb_average(double & grain_size_average);
  void cluster_faces( int start_ilocal, vector<bool> & site_included, int & grain_vol, int & face_sites, vector<int> & faces,
							 bool & multiproc, vector<int> & neigh_procs, vector<double> & cm );
//  void cluster_faces( int start_ilocal, vector<bool> & site_included, int & grain_vol, int & face_sites, vector<int> & faces,
//							 bool & multiproc, vector<int> & neigh_procs );
							 
  void consolidate_mass_center_distributed_grains( vector<int> & spins, vector<double> & cm_list, const int LOCALNUMGRAINS );							 
 
  // Calculating and/or reporting statistics		
  void stats_header(char *strtmp);
  void stats(char *strtmp);

  // Checking lattice configurations
  double count_vacant();
  double count_grain_sites();
  double calculate_density();


 private:
 
  std::map<int,int> hash;
 
  
  vector<int> unique;
  vector<int> adjacent_spins;
  vector<int> grain_start;
  
  vector<int> annihilist_dist;
  vector<int> displacement;
  
  vector<int> annihilation_list;
  int size_annihilist;
  vector<int> annihilation_spin;
  
  vector<bool> check; // Array to check if annihilation is pending
  
  int size_collapsinglist;
  vector<double> collapsing_directions;
	
  vector<int> grain_spins;
  vector<double> grain_mass_center;
  int numgrains;
	
  double Dx, Dy, Dz;
  double xgrid_proc;
  double ygrid_proc;
  double zgrid_proc;
  
  int nx_procs, ny_procs, nz_procs;
  
//  double fraction;

  int frame_depth;
  int border_depth;
  
  int dimension; // Dimension of simulation space
  int nbasis;	// Number of basis atoms in unit cell

  // For density calculations
  int xstart_density, xend_density;
  int ystart_density, yend_density;
  int zstart_density, zend_density;

  // Parameters
  double grain_pore_migration_ratio;
  double pore_migration_annihilation_ratio;
  double pore_migration_temperature; // Temperature for pore migration
  double t_pm_inverse;
  double annihilation_temperature; // Temperature for making vacancy and annihilation
  double t_anni_inverse;
  double time_sinter_start;
  
  double gg_frequency;
  double pm_frequency;
  double a_frequency;
  double gb_factor;
  double gb_ini;
  double density;
  double delta_t;
  	 
  // Info	 
  int vac_made;
  double gsize_ini; // Grain size average at the beginning of sintering
 
  protected:
	int *spin;
};

}

#endif
#endif
