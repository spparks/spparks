/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_MIGRATION_H
#define APP_MIGRATION_H

#include "app_lattice3d.h"

namespace SPPARKS {

class AppMigration : public AppLattice3d {
 public:
  AppMigration(class SPK *, int, char **);
  ~AppMigration();

  double site_energy(int, int, int);
  int site_pick_random(int, int, int, double);		// FOR SWEEPER
  int site_pick_local(int, int, int, double);		// FOR SWEEPER
  double site_propensity(int, int, int, int);
  void site_event(int, int, int, int);
  void site_update_ghosts(int, int, int);
  void site_clear_mask(char ***, int, int, int) {}	// FOR SWEEPER

  void xyz(int, int, int, double *, double *, double *);
  void box_bounds(double *, double *, double *, double *, double *, double *);

 private:
  int nonns;				// number of nearest neighbors (NNs)
  double latt_const;		// lattice constant for the material
  double em_bias;			// applied electromigration bias in +x direction
  double **rates;			// 2D array of hopping rates; size=(nonns+1)x(nonns+1)
							// e.g., rates[4][3]=hopping rate from site w/ 4 NNs to site w/ 3 NNs
  int x1type,x2type;
  int y1type,y2type;		// x,y,z boundary type flags
  int z1type,z2type;
  double vac_bulk;          // concentration of vacancies in bulk and surface/interface layer 
  double vac_surf;          // e.g., num_bulk_atoms*vac_bulk = # empty lattice sites in bulk
							// from 0.0=no vacancies to 1.0=all are vacant 
  int nx1_surf,nx2_surf;
  int ny1_surf,ny2_surf;    // x,y,z surface thicknesses in layers of atoms
  int nz1_surf,nz2_surf;
  int fixed_surf;			// number of layers added to lattice size
  int free_surf;			// for fixed and free surfaces

  // Pointer-to-member functions for different lattice stencils
  void (AppMigration::*nn_fp)(int *, int, int, int);
  double (AppMigration::*bias_fp)(int, int, int, int, int, int);
  void (AppMigration::*xyz_fp)(int,int,int,double *,double *,double *);
  void (AppMigration::*box_fp)(double *,double *,double *);

  // Pointer-to-member functions for different materials
  double (AppMigration::*es_fp)(int);
  void (AppMigration::*cr_fp)();

  void init_app();
  void input_app(char *, int, char **);
  void init_sites();
  void set_boundaries(int, char **);
  void set_vacancies(int, char **);
  void set_bias(int, char **);
  
  int calc_coord_num(int, int, int);
  double hop_propensity(int, int, int, int, int, int);

  // Site event functions for all lattice stencils
  int move_atom(int, int, int, int);
  void broken_test(int, int, int, int, int, int);
  void update_propensity(int, int, int, int, int, int, int);

  // Energy functions for different metals are based on coordination number
  double energy_site_12nn_Al(int);
  void calc_rates_12nn_Al();

  // Functions to get nearest neighbor sites for different lattices
  // first direction is along z, second direction is along x (if an option)
  void find_nnsites_fcc_12nn_100_011(int *, int, int, int);
  void find_nnsites_fcc_12nn_111(int *, int, int, int);
  void find_nnsites_hcp_12nn_0001(int *, int, int, int);

  // Functions to calculate EM force effect on specific neighbor sites
  double calc_bias_fcc_12nn_100_011(int,int,int,int,int,int);
  double calc_bias_fcc_12nn_111(int,int,int,int,int,int);
  double calc_bias_hcp_12nn_0001(int,int,int,int,int,int);

  // Functions to convert a lattice site (i,j,k) to coordinates (x,y,z)
  void get_xyz_fcc_12nn_100_011(int,int,int,double *,double *,double *);
  void get_xyz_fcc_12nn_111(int,int,int,double *,double *,double *);
  void get_xyz_hcp_12nn_0001(int,int,int,double *,double *,double *);
  
  // Functions to calculate upper box bounds for visualization
  void box_bounds_fcc_12nn_100_011(double *, double *, double *);
  void box_bounds_fcc_12nn_111(double *, double *, double *);
  void box_bounds_hcp_12nn_0001(double *, double *, double *);

  // Functions to convert between coordinates and site numbers
  // (for all lattice sites, including ghost points)
  void gsite2ijk(const int&, int&, int&, int&) const;
  int ijk2gsite(const int&, const int&, const int&) const;
};

// Convert from (i,j,k) to site number on my proc; works for ghost points
inline void AppMigration::gsite2ijk(const int& gsite, int& i, int& j, int& k) const {
  i = gsite/((ny_local+2*delghost)*(nz_local+2*delghost)) - delghost+1;
  j = gsite/(nz_local+2*delghost) % (ny_local+2*delghost) - delghost+1;
  k = gsite % (nz_local+2*delghost) - delghost+1;
}

// Convert from site number on my proc to (i,j,k); works for ghost points
inline int AppMigration::ijk2gsite(const int& i, const int& j, const int& k) const {
  return ((i+delghost-1)*((ny_local+2*delghost)*(nz_local+2*delghost)))+((j+delghost-1)*(nz_local+2*delghost))+k+delghost-1;
}
// Arrays for locally owned site numbers are defined in parent class
// NOTE: SITE COORDINATES (i,j,k) ARE UNIVERSAL!!
}

#endif
