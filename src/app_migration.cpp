/* ----------------------------------------------------------------------
SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "app_migration.h"
#include "solve.h"
#include "comm_lattice3d.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace SPPARKS;
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{FIXED,FREE,PERIODIC};
enum{DO_NOT_USE,ATOM,FIXED_ATOM,HOLE,FIXED_HOLE};

/* ---------------------------------------------------------------------- */

AppMigration::AppMigration(SPK *spk, int narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{
	//printf("app_migration running on %d processors\n",nprocs);
	cout << "**********" << endl;
	nonns = 0;
	latt_const = 0.0;
	em_bias = 0.0;
	rates = NULL;
	x1type = -1;
	x2type = -1;
	y1type = -1;
	y2type = -1;
	z1type = -1;
	z2type = -1;
	vac_bulk = 0.0;
	vac_surf = 0.0;
	nx1_surf = 0;
	nx2_surf = 0;
	ny1_surf = 0;
	ny2_surf = 0;
	nz1_surf = 0;
	nz2_surf = 0;
	fixed_surf = 1;
	free_surf = 5;
	nn_fp = NULL;
	es_fp = NULL;
	cr_fp = NULL;
	bias_fp = NULL;
	xyz_fp = NULL;
	box_fp = NULL;

	if (narg != 8) error->all("Invalid app_style migration command");
		
	int iarg=1;
	if (strcmp(arg[iarg],"fcc_12nn_100_011") == 0) {
		nonns = 12;
		nn_fp = &AppMigration::find_nnsites_fcc_12nn_100_011;
		xyz_fp = &AppMigration::get_xyz_fcc_12nn_100_011;
		bias_fp = &AppMigration::calc_bias_fcc_12nn_100_011;
		box_fp = &AppMigration::box_bounds_fcc_12nn_100_011;
	} else if (strcmp(arg[iarg],"fcc_12nn_111") == 0) {
		nonns = 12;
		nn_fp = &AppMigration::find_nnsites_fcc_12nn_111;
		xyz_fp = &AppMigration::get_xyz_fcc_12nn_111;
		bias_fp = &AppMigration::calc_bias_fcc_12nn_111;
		box_fp = &AppMigration::box_bounds_fcc_12nn_111;
	} else if (strcmp(arg[iarg],"hcp_12nn_0001") == 0) {
		nonns = 12;
		nn_fp = &AppMigration::find_nnsites_hcp_12nn_0001;
		xyz_fp = &AppMigration::get_xyz_hcp_12nn_0001;
		bias_fp = &AppMigration::calc_bias_hcp_12nn_0001;
		box_fp = &AppMigration::box_bounds_hcp_12nn_0001;
	} else if (strcmp(arg[iarg],"fcc_12nn_100_010") == 0) {
		error->all("fcc_12nn_100_010 is a work in progress");
	} else if (strcmp(arg[iarg],"fcc_12nn_110") == 0) {
		error->all("fcc_12nn_110 is a work in progress");
	} else error->all("Illegal lattice specifier in app_style migration command");
	cout << "Lattice: " << arg[iarg];
	iarg++;

	if (strcmp(arg[iarg],"Al") == 0) {
		es_fp = &AppMigration::energy_site_12nn_Al;
		cr_fp = &AppMigration::calc_rates_12nn_Al;
	} else error->all("Illegal material specifier in app_style migration command");
	iarg++;
	
	latt_const = atof(arg[iarg++]);
	nx_global = atoi(arg[iarg++]);
	ny_global = atoi(arg[iarg++]);
	nz_global = atoi(arg[iarg++]);
	cout << "; Nx = " << nx_global << ", Ny = " << ny_global << ", Nz = " << nz_global << endl;
	//printf("Lattice constant = %g\n",latt_const);
	
	seed = atoi(arg[iarg++]);
	random = new RandomPark(seed);
	
	dellocal = 1;
	delghost = 2;

	// Create the hopping rates matrix
	if (rates == NULL)
		memory->create_2d_T_array(rates,(nonns+1),(nonns+1),"app_migration:rates");
	
}

/* ---------------------------------------------------------------------- */

AppMigration::~AppMigration()
{
	delete random;
    memory->destroy_3d_T_array(lattice,nxlo,nylo,nzlo);
	delete comm;
	memory->destroy_2d_T_array(rates);
}

/* ---------------------------------------------------------------------- */

void AppMigration::init_app()
{
	(this->*cr_fp)();		// calculate hopping rates

	// Modify lattice size to account for certain boundary conditions
	// Fixed boundaries => add fixed_surf atom layers to lattice size
	// Free boundaries => add free_surf atom layers to lattice size
	
	if (x1type == FIXED) nx_global += fixed_surf;
	else if (x1type == FREE) nx_global += free_surf;
	else if (x1type == -1) error->all("x1type not set yet");
	
	if (x2type == FIXED) nx_global += fixed_surf;
	else if (x2type == FREE) nx_global += free_surf;
	else if (x2type == -1) error->all("x2type not set yet");
	
	if (y1type == FIXED) ny_global += fixed_surf;
	else if (y1type == FREE) ny_global += free_surf;
	else if (y1type == -1) error->all("y1type not set yet");
	
	if (y2type == FIXED) ny_global += fixed_surf;
	else if (y2type == FREE) ny_global += free_surf;
	else if (y2type == -1) error->all("y2type not set yet");
	
	if (z1type == FIXED) nz_global += fixed_surf;
	else if (z1type == FREE) nz_global += free_surf;
	else if (z1type == -1) error->all("z1type not set yet");
	
	if (z2type == FIXED) nz_global += fixed_surf;
	else if (z2type == FREE) nz_global += free_surf;
	else if (z2type == -1) error->all("z2type not set yet");
	
	// define lattice and partition it across processors
	procs2lattice();
	memory->create_3d_T_array(lattice,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
				  "appmigration:lattice");
	
	// initialize my portion of lattice sites
	init_sites();

	// Output for testing
	cout << "Boundary types are (0=fixed, 1=free, 2=periodic): " << x1type << x2type << y1type << y2type << z1type << z2type << endl;
	cout << "nglobal for (x,y,z) on proc "<< me << " are " << nx_global <<","<<ny_global<<","<<nz_global<< endl;			
	cout << "offsets for (x,y,z) on proc "<< me << " are " << nx_offset <<","<<ny_offset<<","<<nz_offset<< endl;
	cout << " nlocal for (x,y,z) on proc "<< me << " are " << nx_local <<","<<ny_local<<","<<nz_local<< endl;			
}

/* ---------------------------------------------------------------------- */

void AppMigration::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"boundaries") == 0) set_boundaries(narg,arg);
  else if (strcmp(command,"vacancies") == 0) set_vacancies(narg,arg);
  else if (strcmp(command,"bias") == 0) set_bias(narg,arg);
  else error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   Initialize free boundary, bulk, and surface vacancies (globally)
------------------------------------------------------------------------- */
void AppMigration::init_sites()
{
	int i,j,k;
	
	// initialize local sites to occupied (ATOM)
	for (i = 1; i <= nx_local; i++)
	for (j = 1; j <= ny_local; j++)
	for (k = 1; k <= nz_local; k++)
		lattice[i][j][k] = ATOM;
	
	// initialize free boundary layers to unoccupied (HOLE)
	if (x1type == FREE && nx_offset == 0)
	{	//proc with x1 global edge
		for (i = 1; i <= free_surf; i++)
		for (j = 1; j <= ny_local; j++)
		for (k = 1; k <= nz_local; k++) 
			lattice[i][j][k] = HOLE;
	}
	if (x2type == FREE && (nx_local+nx_offset == nx_global) )
	{	//proc with x2 global edge
		for (i = nx_local-free_surf+1; i <= nx_local; i++)
		for (j = 1; j <= ny_local; j++)
		for (k = 1; k <= nz_local; k++)
			lattice[i][j][k] = HOLE;
	}
	if (y1type == FREE && ny_offset == 0)
	{	//proc with y1 global edge
		for (i = 1; i <= nx_local; i++)
		for (j = 1; j <= free_surf; j++)
		for (k = 1; k <= nz_local; k++)
			lattice[i][j][k] = HOLE;
	}
	if (y2type == FREE && (ny_local+ny_offset == ny_global) )
	{	//proc with y2 global edge
		for (i = 1; i <= nx_local; i++)
		for (j = ny_local-free_surf+1; j <= ny_local; j++)		
		for (k = 1; k <= nz_local; k++)
			lattice[i][j][k] = HOLE;
	}
	if (z1type == FREE && nz_offset == 0)
	{	//proc with z1 global edge
		for (i = 1; i <= nx_local; i++)
		for (j = 1; j <= ny_local; j++) 
		for (k = 1; k <= free_surf; k++)
			lattice[i][j][k] = HOLE;
	}
	if (z2type == FREE && (nz_local+nz_offset == nz_global) )
	{	//proc with z2 global edge
		for (i = 1; i <= nx_local; i++)
		for (j = 1; j <= ny_local; j++) 
		for (k = nz_local-free_surf+1; k <= nz_local; k++)
			lattice[i][j][k] = HOLE;
	}

	// set global edge layers to FIXED_HOLE for free 
	// or FIXED_ATOM for fixed boundaries
	if (x1type == FIXED || x1type == FREE) {
		if (nx_offset == 0) {
			i = 1;
			for (j = 1; j <= ny_local; j++) {
			for (k = 1; k <= nz_local; k++) {
				if (x1type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (x1type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}
	if (x2type == FIXED || x2type == FREE) {
		if (nx_local+nx_offset == nx_global) {
			i = nx_local;
			for (j = 1; j <= ny_local; j++) {
			for (k = 1; k <= nz_local; k++) {
				if (x2type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (x2type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}
	if (y1type == FIXED || y1type == FREE) {
		if (ny_offset == 0) {
			j = 1;
			for (i = 1; i <= nx_local; i++) {
			for (k = 1; k <= nz_local; k++) {
				if (y1type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (y1type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}
	if (y2type == FIXED || y2type == FREE) {
		if (ny_local+ny_offset == ny_global) {
			j = ny_local;
			for (i = 1; i <= nx_local; i++) {
			for (k = 1; k <= nz_local; k++) {
				if (y2type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (y2type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}
	if (z1type == FIXED || z1type == FREE) {
		if (nz_offset == 0) {
			k = 1;
			for (i = 1; i <= nx_local; i++) {
			for (j = 1; j <= ny_local; j++) {
				if (z1type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (z1type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}
	if (z2type == FIXED || z2type == FREE) {
		if (nz_local+nz_offset == nz_global) {
			k = nz_local;
			for (i = 1; i <= nx_local; i++) {
			for (j = 1; j <= ny_local; j++) {
				if (z2type == FIXED) lattice[i][j][k] = FIXED_ATOM;
				else if (z2type == FREE) lattice[i][j][k] = FIXED_HOLE;
	}}}}

	// determine interior coordinates of global lattice
	// (not counting additional fixed or free surfaces)
	// for purpose of assigning vacancies
	int x1,x2,y1,y2,z1,z2;

	if (x1type == FREE) x1 = free_surf+1;
	else if (x1type == FIXED) x1 = fixed_surf+1;
	else x1 = 1;
	
	if (x2type == FREE) x2 = nx_global - free_surf;
	else if (x2type == FIXED) x2 = nx_global - fixed_surf;
	else x2 = nx_global;
	
	if (y1type == FREE) y1 = free_surf+1;
	else if (y1type == FIXED) y1 = fixed_surf+1;
	else y1 = 1;
	
	if (y2type == FREE) y2 = ny_global - free_surf;
	else if (y2type == FIXED) y2 = ny_global - fixed_surf;
	else y2 = ny_global;

	if (z1type == FREE) z1 = free_surf+1;
	else if (z1type == FIXED) z1 = fixed_surf+1;
	else z1 = 1;
	
	if (z2type == FREE) z2 = nz_global - free_surf;
	else if (z2type == FIXED) z2 = nz_global - fixed_surf;
	else z2 = nz_global;
	
	// calculate number of GLOBAL bulk and surface lattice sites
	int n_bulk = (x2-x1-nx1_surf-nx2_surf)*(y2-y1-ny1_surf-ny2_surf)*(z2-z1-nz1_surf-nz2_surf);
	int n_surf = ((x2-x1)*(y2-y1)*(z2-z1))-n_bulk;

	// calculate number of LOCAL bulk and surface vacancies
	int n_bulk_vac = (int)floor(vac_bulk*((double)n_bulk)/nprocs);
	int n_surf_vac = (int)floor(vac_surf*((double)n_surf)/nprocs);

	// loop over global list so that assigment is independent 
	// of parallel decomposition and also so that each local 
	// domain is initialized with different vacant sites.

	// initialize bulk vacancies 
	int ii,jj,kk;
	if (n_bulk_vac > 0) {
		int bulk_counter = 0;
		// calculate bulk boundaries on my proc
		int i1 = MAX(x1 + nx1_surf, 1 + nx_offset);
		int i2 = MIN(x2 - nx2_surf, nx_local + nx_offset);
		int j1 = MAX(y1 + ny1_surf, 1 + ny_offset);
		int j2 = MIN(y2 - ny2_surf, ny_local + ny_offset);
		int k1 = MAX(z1 + nz1_surf, 1 + nz_offset);
		int k2 = MIN(z2 - nz2_surf, nz_local + nz_offset);
		
		while (bulk_counter < n_bulk_vac) {
			// pick random indices to give a global bulk site on my proc
			i = i1 + random->irandom(i2-i1+1)-1;
			j = j1 + random->irandom(j2-j1+1)-1;
			k = k1 + random->irandom(k2-k1+1)-1;
			// convert global to local indices
			ii = i - nx_offset;
			jj = j - ny_offset;
			kk = k - nz_offset;
			// make site vacant if it isn't already
			if (lattice[ii][jj][kk] != HOLE) {
				lattice[ii][jj][kk] = HOLE;
				bulk_counter++;
			}
		}
	}

	// need to assign surface vacancies (see deadcode.cpp for previous code)

	// equilibrate system for ??? seconds before turning on bias

}

/* ---------------------------------------------------------------------- */

void AppMigration::set_boundaries(int narg, char **arg)
{
	if (narg != 6) error->all("Illegal boundaries command");

	if (strcmp(arg[0],"fixed") == 0) x1type = FIXED;
	else if (strcmp(arg[0],"free") == 0) x1type = FREE;
	else if (strcmp(arg[0],"periodic") == 0) x1type = PERIODIC;
	else error->all("Invalid x1 boundary type entered");
	
	if (strcmp(arg[1],"fixed") == 0) x2type = FIXED;
	else if (strcmp(arg[1],"free") == 0) x2type = FREE;
	else if (strcmp(arg[1],"periodic") == 0) x2type = PERIODIC;
	else error->all("Invalid x2 boundary type entered");

	if (strcmp(arg[2],"fixed") == 0) y1type = FIXED;
	else if (strcmp(arg[2],"free") == 0) y1type = FREE;
	else if (strcmp(arg[2],"periodic") == 0) y1type = PERIODIC;
	else error->all("Invalid y1 boundary type entered");

	if (strcmp(arg[3],"fixed") == 0) y2type = FIXED;
	else if (strcmp(arg[3],"free") == 0) y2type = FREE;
	else if (strcmp(arg[3],"periodic") == 0) y2type = PERIODIC;
	else error->all("Invalid y2 boundary type entered");

	if (strcmp(arg[4],"fixed") == 0) z1type = FIXED;
	else if (strcmp(arg[4],"free") == 0) z1type = FREE;
	else if (strcmp(arg[4],"periodic") == 0) z1type = PERIODIC;
	else error->all("Invalid z1 boundary type entered");

	if (strcmp(arg[5],"fixed") == 0) z2type = FIXED;
	else if (strcmp(arg[5],"free") == 0) z2type = FREE;
	else if (strcmp(arg[5],"periodic") == 0) z2type = PERIODIC;
	else error->all("Invalid z2 boundary type entered");
}

/* ---------------------------------------------------------------------- */

void AppMigration::set_vacancies(int narg, char **arg)
{
	if (narg != 8) error->all("Illegal vacancies command");
	vac_bulk = atof(arg[0]);
	vac_surf = atof(arg[1]);
	nx1_surf = atoi(arg[2]);
	nx2_surf = atoi(arg[3]);
	ny1_surf = atoi(arg[4]);
	ny2_surf = atoi(arg[5]);
	nz1_surf = atoi(arg[6]);
	nz2_surf = atoi(arg[7]);
}

/* ---------------------------------------------------------------------- */

void AppMigration::set_bias(int narg, char **arg)
{
	if (narg != 1) error->all("Illegal bias command");
	em_bias = atof(arg[0]);
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */
double AppMigration::site_energy(int i, int j, int k)
{
	int cn = 0;
	double eng = 0.0;

	// if unoccupied site, energy is 0; if occupied site, calculate energy
	
	if (lattice[i][j][k] == ATOM) {					// mobile atom
		cn = calc_coord_num(i,j,k);
		eng = (this->*es_fp)(cn);
	} else if (lattice[i][j][k] == FIXED_ATOM) {	// global boundary atom
		int ii,jj,kk,neighbors[nonns];		
		(this->*nn_fp)(neighbors,i,j,k);
		for (int n = 0; n < nonns; n++) {
			gsite2ijk(neighbors[n],ii,jj,kk);
			if (ii < 1 || ii > nx_global ||			// treat non-local "neighbors"
				jj < 1 || jj > ny_global ||			// of fixed global boundary 
				kk < 1 || kk > nz_global ) cn++;	// sites as occupied
			else	// check occupation for local neighbors
				if (lattice[ii][jj][kk] == ATOM || lattice[ii][jj][kk] == FIXED_ATOM)
					cn++;				
		}
		eng = (this->*es_fp)(cn);
	}

	return eng;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly - FOR SWEEPER
------------------------------------------------------------------------- */
int AppMigration::site_pick_random(int i, int j, int k, double ran)
{
  int iran = 1;
  return iran;
}

/* ----------------------------------------------------------------------
   pick new state for site randomly from neighbor values - FOR SWEEPER
------------------------------------------------------------------------- */
int AppMigration::site_pick_local(int i, int j, int k, double ran)
{
  int iran = 1;
  return iran;
}

/* ----------------------------------------------------------------------
   compute total propensity of local site
------------------------------------------------------------------------- */
double AppMigration::site_propensity(int i, int j, int k, int full)
{
	// If proc owns full domain and we have PBCs, update ghost sites
	if (full) {
		if (x1type == PERIODIC || y1type == PERIODIC || z1type == PERIODIC)
			site_update_ghosts(i,j,k);
	}

	// Total propensity for site (i,j,k) is a sum of the individual
	// propensities for a hop from (i,j,k) to each of its neighbors
	double prop = 0.0;
	int n,ii,jj,kk;
	int neighbors[nonns];
	
	(this->*nn_fp)(neighbors,i,j,k);
	for (n = 0; n < nonns; n++) {
		gsite2ijk(neighbors[n],ii,jj,kk);
		prop += hop_propensity(i,j,k,ii,jj,kk);
	}
	return prop;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, neighbor sites may be across PBC
   (Note: ghost site updates take place in called functions)
   if proc owns sector, ignore neighbor sites outside sector
------------------------------------------------------------------------- */
void AppMigration::site_event(int i, int j, int k, int full)
{		
	// Event: 1) occupied site (i,j,k) chosen by solver becomes unoccupied
	//		  2) empty neighbor site is selected to become occupied based 
	//			 on individual site propensities
	int gsite = move_atom(i,j,k,full);
	int ii,jj,kk;
	gsite2ijk(gsite,ii,jj,kk);

	// Test to see if the wire is broken - DUMMY FUNCTION FOR NOW
	broken_test(i,j,k,ii,jj,kk);

	// update propensities for (i,j,k), (ii,jj,kk), and affected neighbors
	update_propensity(i,j,k,ii,jj,kk,full);

}

/* ----------------------------------------------------------------------
   update neighbors of site if neighbors are ghost cells -- called by 
   move_atom() and site_propensity() when single proc owns entire domain
------------------------------------------------------------------------- */
void AppMigration::site_update_ghosts(int i, int j, int k)
{
	int ii,jj,kk;
	
	// NOTE: periodic boundary conditions must occur in pairs to make sense.
	// Thus, if x1type is PERIODIC, x2type must be also, and same for y and z.
	
	if (x1type == PERIODIC) {
		if (i == 1) {
			for (jj = j-1; jj <= j+1; jj++)
				for (kk = k-1; kk <= k+1; kk++)
					lattice[i-1][jj][kk] = lattice[nx_local][jj][kk];
		}
		if (i == nx_local) {
			for (jj = j-1; jj <= j+1; jj++)
				for (kk = k-1; kk <= k+1; kk++)
					lattice[i+1][jj][kk] = lattice[1][jj][kk];
		}
	}

	if (y1type == PERIODIC) {
		if (j == 1) {
			for (ii = i-1; ii <= i+1; ii++)
				for (kk = k-1; kk <= k+1; kk++)
					lattice[ii][j-1][kk] = lattice[ii][ny_local][kk];
		}
		if (j == ny_local) {
			for (ii = i-1; ii <= i+1; ii++)
				for (kk = k-1; kk <= k+1; kk++)
					lattice[ii][j+1][kk] = lattice[ii][1][kk];
		}
	}

	if (z1type == PERIODIC) {
		if (k == 1) {
			for (ii = i-1; ii <= i+1; ii++)
				for (jj = j-1; jj <= j+1; jj++)
					lattice[ii][jj][k-1] = lattice[ii][jj][nz_local];
		}
		if (k == nz_local) {
			for (ii = i-1; ii <= i+1; ii++)
				for (jj = j-1; jj <= j+1; jj++)
					lattice[ii][jj][k+1] = lattice[ii][jj][1];
		}
	}
	
}

/* ----------------------------------------------------------------------
   Calculate the coordination number of an atom at site (i,j,k)
------------------------------------------------------------------------- */
int AppMigration::calc_coord_num(int i, int j, int k)
{
	int nocc = 0;
	
	int n,ii,jj,kk,neighbors[nonns];
	(this->*nn_fp)(neighbors,i,j,k);			// get nearest neighbor sites
	for (n = 0; n < nonns; n++) {
		gsite2ijk(neighbors[n],ii,jj,kk);		// convert site to coordinates
		if (lattice[ii][jj][kk] == ATOM ||		// count occupied neighbors
			lattice[ii][jj][kk] == FIXED_ATOM) 
			nocc++;
	}
	return nocc;
}

/* ----------------------------------------------------------------------
   Calculate the propensity for a hop from (i,j,k) to (ii,jj,kk)
   ~ propensity = 0 for a site (i,j,k) that is fixed and/or unoccupied
					and for a hop-to site (ii,jj,kk) that is not unoccupied
   ~ propensity is calculated for an occupied, non-fixed site (ATOM) (i,j,k) 
     that can hop to an unoccupied, non-fixed site (HOLE) (ii,jj,kk)
   ~ fixed sites relate to the global boundary conditions:
	   fixed_surf layers of sites on fixed surfaces (sites that remain occupied)
	   last layer of sites on free surfaces (so atoms can't hop off lattice)
------------------------------------------------------------------------- */
double AppMigration::hop_propensity(int i, int j, int k, int ii, int jj, int kk)
{
	double hopprop;

	if (lattice[i][j][k] == ATOM) {			// (i,j,k) is occupied & not fixed
		if (lattice[ii][jj][kk] == HOLE) {	// (ii,jj,kk) is unoccupied & not fixed
			// move is valid ==> determine hopping propensity
			int cn_init = calc_coord_num(i,j,k);
			int cn_final = calc_coord_num(ii,jj,kk)-1;
			if (cn_final < 0) cn_final = 0;
			hopprop = (this->*bias_fp)(i,j,k,ii,jj,kk)*rates[cn_init][cn_final];
		} else hopprop = 0.0;
	} else hopprop = 0.0;
	
	return hopprop;
}

/* ----------------------------------------------------------------------
   Check area of atom just moved to see if wire is broken
------------------------------------------------------------------------- */
void AppMigration::broken_test(int i, int j, int k, int ii, int jj, int kk) 
{
}

/* ----------------------------------------------------------------------
   Move atom chosen by the solver to an empty neighboring site
------------------------------------------------------------------------- */
int AppMigration::move_atom(int i, int j, int k, int full)
{
	// If proc owns full domain and we have PBCs, update ghost sites
	if (full) {
		if (x1type == PERIODIC || y1type == PERIODIC || z1type == PERIODIC)
			site_update_ghosts(i,j,k);
	}

	int ii,jj,kk,n,ncand,candlist[nonns],neighbors[nonns];
	double testprop,candprop[nonns];
	double totalprop = 0.0;
	
	// Find empty nearest neighbor sites (candidate sites for the move)
	(this->*nn_fp)(neighbors,i,j,k);			// get nearest neighbor sites
	ncand = 0;
	for (n = 0; n < nonns; n++) {
		gsite2ijk(neighbors[n],ii,jj,kk);		// convert site to coordinates
		testprop = hop_propensity(i,j,k,ii,jj,kk);	// calc hopping propensity
		if (testprop > 0.0) {					// if move is possible (prop>0)
			candlist[ncand] = neighbors[n];		// add site to candidate list,
			candprop[ncand] = testprop;			// store propensity, add it to
			totalprop += candprop[ncand];		// the total propensity, and 
			ncand++;							// increment candidate counter
		}
	}
	
	//***Error output if no candidate sites are found
	if (ncand == 0) {
		printf("At timestep %d (time = %g), atom to move is (%d,%d,%d), site %d (gsite %d)\n",
				ntimestep,time,i,j,k,ijk2site[i][j][k],ijk2gsite(i,j,k));
		printf("atom has propensity = %g, site_prop = %g, totalprop = %g\n",
				propensity[ijk2site[i][j][k]],site_propensity(i,j,k,full),totalprop);
		error->all("Candidate site not found!");
	} //***end error output

	// Pick site to move to randomly based on propensity
	double randnum = totalprop * random->uniform();
	double sum = 0.0;
	for (n = 0; n < ncand; n++) {
		sum += candprop[n];
		if (sum > randnum) break;
	}

	//***Error output if none of the candidate sites are chosen
	if (n >= ncand) {
		printf("At timestep %d (time = %g), atom to move is (%d,%d,%d), site %d (gsite %d)\n",
				ntimestep,time,i,j,k,ijk2site[i][j][k],ijk2gsite(i,j,k));
		printf("atom has propensity = %g, site_prop = %g, totalprop = %g\n",
				propensity[ijk2site[i][j][k]],site_propensity(i,j,k,full),totalprop);
		printf("Chosen n = %d >= number of candidates = %d! sum = %g, randnum = %g\n",
				n,ncand,sum,randnum);
		for (int m = 0; m < ncand; m++) {
			printf("candlist[%d] = %d, candprop[%d] = %g\n",m,candlist[m],m,candprop[m]);
		}
		error->all("No candidate site chosen!");
	} //***end error output

	// Make the move
	gsite2ijk(candlist[n],ii,jj,kk);	// get coordinates of site to jump to
	lattice[ii][jj][kk] = ATOM;			// site (ii,jj,kk) becomes occupied
	lattice[i][j][k] = HOLE;			// site (i,j,k) becomes unoccupied

	return ijk2gsite(ii,jj,kk);
}

/* ----------------------------------------------------------------------
   Update propensities for all affected sites (locally owned)
------------------------------------------------------------------------- */
void AppMigration::update_propensity(int i,int j,int k,int ii,int jj,int kk,int full)
{	
	// Check if we have any PBCs
	bool pbc;
	if (x1type == PERIODIC || y1type == PERIODIC || z1type == PERIODIC)
		pbc = true;
	else pbc = false;

	// Sites that need a propensity update include the neighbors of the
	// neighbors of sites (i,j,k) and (ii,jj,kk). Since every lattice 
	// stencil contains the nearest neighbors within +/- 1 of i,j,k, 
	// it is sufficient to check only lattice sites within +/- 2 of
	// i,j,k and ii,jj,kk for propensity changes.
	int iii,jjj,kkk,iloop,jloop,kloop,flag,tempsite,temp_sites[150];
	int count = 0;			// count number of sites to update
	double tempprop;

	for (iloop = MIN(i,ii)-2; iloop <= MAX(i,ii)+2; iloop++)
	for (jloop = MIN(j,jj)-2; jloop <= MAX(j,jj)+2; jloop++)
	for (kloop = MIN(k,kk)-2; kloop <= MAX(k,kk)+2; kloop++) {
		iii = iloop; jjj = jloop; kkk = kloop;
		flag = 1;
		if (full && pbc) ijkpbc(iii,jjj,kkk);
		 // update for my sector's sites only
		else if ( iii < nx_sector_lo || iii > nx_sector_hi ||
				  jjj < ny_sector_lo || jjj > ny_sector_hi || 
				  kkk < nz_sector_lo || kkk > nz_sector_hi ) flag = 0;
		if (flag) {	 
			tempsite = ijk2site[iii][jjj][kkk];
			tempprop = site_propensity(iii,jjj,kkk,full);
			if (propensity[tempsite] != tempprop) {	// site needs update
				temp_sites[count] = tempsite;
				propensity[tempsite] = tempprop;
				count++;
			}
		}
	}

	// Create array of sites to be updated to pass to solver
	if (count > 0) {
		int update_sites[count];
		for (int n = 0; n < count; n++) {
			update_sites[n] = temp_sites[n];
		}
		solve->update(count,update_sites,propensity);
	}
}

/* ---------------------------------------------------------------------- */
// MATERIAL-SPECIFIC FUNCTIONS
/* ----------------------------------------------------------------------
   Calculate energy of a site as a function of coordination number for Al
------------------------------------------------------------------------- */
double AppMigration::energy_site_12nn_Al(int cn)
{
	int n;
	double f[13], E[13];
	double phi = -0.15;

	f[0] =  0.0;
	f[1] =  0.0;
	f[2] =  0.0;
	f[3] = -1.960;
	f[4] = -2.044;
	f[5] = -2.153;
	f[6] = -2.216;
	f[7] = -2.248;
	f[8] = -2.284;
	f[9] = -2.303;
	f[10]= -2.305;
	f[11]= -2.378;
	f[12]= -2.460;

	for (n = 0; n < 13; n++) E[n] = f[n] + 0.5*phi*(double)n;
	
	return E[cn];		// in eV
}

/* ----------------------------------------------------------------------
   Compute the hopping rate matrix for jumps from site i to site j in Al
------------------------------------------------------------------------- */
void AppMigration::calc_rates_12nn_Al()
{
	int i,j;
	double Ei, Ej;
	double Emig[13], prefactor[13], Ediff[13][13];
	double kT = 0.000086173*temperature;	// in eV
	
	Emig[0] = 0.0;		// in eV
	Emig[1] = 0.0;
	Emig[2] = 0.0;
	Emig[3] = 0.10;
	Emig[4] = 0.33;
	Emig[5] = 0.33;
	Emig[6] = 0.60;
	Emig[7] = 0.60;
	Emig[8] = 0.60;
	Emig[9] = 0.60;
	Emig[10]= 0.60;
	Emig[11]= 0.60;
	Emig[12]= 0.0;
	
	prefactor[0] = 0.0;		// in cm^2 / s
	prefactor[1] = 0.0;
	prefactor[2] = 0.0;
	prefactor[3] = 0.0002;
	prefactor[4] = 0.4;
	prefactor[5] = 0.004;
	prefactor[6] = 0.01;
	prefactor[7] = 0.01;
	prefactor[8] = 0.01;
	prefactor[9] = 0.01;
	prefactor[10]= 0.01;
	prefactor[11]= 0.01;
	prefactor[12]= 0.0;
	
	for (i = 0; i < 13; i++) {
		for (j = 0; j < 13; j++) {
			if (j < 3) rates[i][j] = 0.0; //can't hop to a site with < 3 NNs
			else {
				Ei = (this->*es_fp)(i);
				Ej = (this->*es_fp)(j);
				if ((Ej - Ei) > 0.0)
					Ediff[i][j] = Ej - Ei;
				else
					Ediff[i][j] = 0.0;
				rates[i][j] = prefactor[i]*exp(-1.0*((Emig[i]+Ediff[i][j])/kT));
			}
		}
	}
}

/* ---------------------------------------------------------------------- */
// LATTICE-SPECIFIC FUNCTIONS
/* ----------------------------------------------------------------------
   Get nearest neighbor sites for site (i,j,k) on fcc(100)[011]
   (100) surface normal to z axis, (011) surface normal to x axis
------------------------------------------------------------------------- */
void AppMigration::find_nnsites_fcc_12nn_100_011(int *nnsites, int i, int j, int k)
{
	int n, ni[12], nj[12], nk[12];

	// Hard code indices for nearest neighbors
	if (k % 2 == 0)		// even layer
	{
		ni[0] = i-1; nj[0] = j-1; nk[0] = k-1;
		ni[1] = i-1; nj[1] = j-1; nk[1] = k+1;
		ni[2] = i-1; nj[2] = j  ; nk[2] = k-1;
		ni[3] = i-1; nj[3] = j  ; nk[3] = k  ;
		ni[4] = i-1; nj[4] = j  ; nk[4] = k+1;
		ni[5] = i  ; nj[5] = j-1; nk[5] = k-1;
		ni[6] = i  ; nj[6] = j-1; nk[6] = k  ;
		ni[7] = i  ; nj[7] = j-1; nk[7] = k+1;
		ni[8] = i  ; nj[8] = j  ; nk[8] = k-1;
		ni[9] = i  ; nj[9] = j  ; nk[9] = k+1;
		ni[10]= i  ; nj[10]= j+1; nk[10]= k  ;
		ni[11]= i+1; nj[11]= j  ; nk[11]= k  ;
	}
	else				// odd layer
	{
		ni[0] = i-1; nj[0] = j  ; nk[0] = k  ;
		ni[1] = i  ; nj[1] = j-1; nk[1] = k  ;
		ni[2] = i  ; nj[2] = j  ; nk[2] = k-1;
		ni[3] = i  ; nj[3] = j  ; nk[3] = k+1;
		ni[4] = i  ; nj[4] = j+1; nk[4] = k-1;
		ni[5] = i  ; nj[5] = j+1; nk[5] = k  ;
		ni[6] = i  ; nj[6] = j+1; nk[6] = k+1;
		ni[7] = i+1; nj[7] = j  ; nk[7] = k-1;
		ni[8] = i+1; nj[8] = j  ; nk[8] = k  ;
		ni[9] = i+1; nj[9] = j  ; nk[9] = k+1;
		ni[10]= i+1; nj[10]= j+1; nk[10]= k-1;
		ni[11]= i+1; nj[11]= j+1; nk[11]= k+1;	
	}

	// Convert neighbor indices to sites
	for (n = 0; n < 12; n++) {
		nnsites[n] = ijk2gsite(ni[n],nj[n],nk[n]);
	}
}

/* ----------------------------------------------------------------------
   Calculate EM bias on nearest neighbor sites of (i,j,k) on fcc(100)[011]
------------------------------------------------------------------------- */
double AppMigration::calc_bias_fcc_12nn_100_011(int i,int j,int k,int ii,int jj,int kk)
{
	double kT = 0.000086173*temperature;
	double delta = em_bias;		// NEEDS TO BE FIGURED OUT
								// global variable em_bias is input value
	double E = 0.0;
	
	// Hard code EM bias for nearest neighbors

	if (kk == k)	// In-plane neighbors (bias is same for all layers)
	{
		if (ii == i+1) E = delta;
		else if (ii == i-1) E = (-1.0)*delta;
	}
	else		// Out-of-plane neighbors (bias depends on which layer)
	{
		if (k % 2 == 0)			// even layer
		{
			if (ii == i) E = (0.5)*delta;
			else if (ii == i-1) E = (-0.5)*delta;
		}
		else					// odd layer
		{
			if (ii == i+1) E = (0.5)*delta;
			else if (ii == i) E = (-0.5)*delta;
		}
	}

	return exp(E/kT);
}

/* ----------------------------------------------------------------------
   Convert site (i,j,k) to physical coordinates (x,y,z) on fcc(100)[011]
------------------------------------------------------------------------- */
void AppMigration::get_xyz_fcc_12nn_100_011(int i, int j, int k, 
											double *x, double *y, double *z)
{
	// global variable latt_const is the material's lattice constant, a

	double sq = 1.0 / sqrt(2.0);

	// global array positions
	int gi = i + nx_offset;
	int gj = j + ny_offset;
	int gk = k + nz_offset;
	
	// x and y values depend on even or odd layer (gk)
	if (gk % 2 == 0)		// even layer
	{
		*x = (double)(gi)*sq*latt_const;
		*y = (double)(gj)*sq*latt_const;
	}
	else					// odd layer
	{
		*x = ((double)(gi)+0.5)*sq*latt_const;
		*y = ((double)(gj)+0.5)*sq*latt_const;
	}

	// z value is the same for all gi, gj
	*z = (double)(gk)*latt_const/2.0;	
}

/* ----------------------------------------------------------------------
   Calculate box upper bounds for visualization on fcc(100)[011]
------------------------------------------------------------------------- */
void AppMigration::box_bounds_fcc_12nn_100_011(double *xhi,double *yhi,double *zhi)
{
	// global variable latt_const is the material's lattice constant, a

	double sq = 1.0 / sqrt(2.0);

	*xhi = ((double)(nx_global) + 0.5)*sq*latt_const;
	*yhi = ((double)(ny_global) + 0.5)*sq*latt_const;
	*zhi = (double)(nz_global)*0.5*latt_const;
}

/* ----------------------------------------------------------------------
   Get nearest neighbor sites for site (i,j,k) on fcc(111)
   (111) surface normal to z axis, (-110) along x axis
------------------------------------------------------------------------- */
void AppMigration::find_nnsites_fcc_12nn_111(int *nnsites, int i, int j, int k)
{
	int n, ni[12], nj[12], nk[12];

	// Hard code indices for nearest neighbors
	if (k % 3 == 0)			// A layer
	{
		ni[0] = i-1; nj[0] = j  ; nk[0] = k  ;
		ni[1] = i-1; nj[1] = j  ; nk[1] = k+1;
		ni[2] = i  ; nj[2] = j-1; nk[2] = k-1;
		ni[3] = i  ; nj[3] = j-1; nk[3] = k  ;
		ni[4] = i  ; nj[4] = j-1; nk[4] = k+1;
		ni[5] = i  ; nj[5] = j  ; nk[5] = k-1;
		ni[6] = i  ; nj[6] = j  ; nk[6] = k+1;
		ni[7] = i  ; nj[7] = j+1; nk[7] = k  ;
		ni[8] = i+1; nj[8] = j-1; nk[8] = k-1;
		ni[9] = i+1; nj[9] = j-1; nk[9] = k  ;
		ni[10]= i+1; nj[10]= j  ; nk[10]= k  ;
		ni[11]= i+1; nj[11]= j+1; nk[11]= k  ;
	}
	else if (k % 3 == 1)	// B layer
	{
		ni[0] = i-1; nj[0] = j  ; nk[0] = k  ;
		ni[1] = i  ; nj[1] = j-1; nk[1] = k  ;
		ni[2] = i  ; nj[2] = j  ; nk[2] = k-1;
		ni[3] = i  ; nj[3] = j  ; nk[3] = k+1;
		ni[4] = i  ; nj[4] = j+1; nk[4] = k  ;
		ni[5] = i+1; nj[5] = j-1; nk[5] = k  ;
		ni[6] = i+1; nj[6] = j-1; nk[6] = k+1;
		ni[7] = i+1; nj[7] = j  ; nk[7] = k-1;
		ni[8] = i+1; nj[8] = j  ; nk[8] = k  ;
		ni[9] = i+1; nj[9] = j  ; nk[9] = k+1;
		ni[10]= i+1; nj[10]= j+1; nk[10]= k-1;
		ni[11]= i+1; nj[11]= j+1; nk[11]= k  ;	
	}
	else if (k % 3 == 2)	// C layer
	{
		ni[0] = i-1; nj[0] = j  ; nk[0] = k  ;
		ni[1] = i-1; nj[1] = j  ; nk[1] = k-1;
		ni[2] = i  ; nj[2] = j-1; nk[2] = k  ;
		ni[3] = i  ; nj[3] = j  ; nk[3] = k-1;
		ni[4] = i  ; nj[4] = j  ; nk[4] = k+1;
		ni[5] = i  ; nj[5] = j+1; nk[5] = k-1;
		ni[6] = i  ; nj[6] = j+1; nk[6] = k  ;
		ni[7] = i  ; nj[7] = j+1; nk[7] = k+1;
		ni[8] = i+1; nj[8] = j-1; nk[8] = k  ;
		ni[9] = i+1; nj[9] = j  ; nk[9] = k  ;
		ni[10]= i+1; nj[10]= j+1; nk[10]= k  ;
		ni[11]= i+1; nj[11]= j+1; nk[11]= k+1;
	}

	// Convert neighbor indices to sites
	for (n = 0; n < 12; n++) {
		nnsites[n] = ijk2gsite(ni[n],nj[n],nk[n]);
	}
}

/* ----------------------------------------------------------------------
   Calculate EM bias on nearest neighbor sites of (i,j,k) on fcc(111)
------------------------------------------------------------------------- */
double AppMigration::calc_bias_fcc_12nn_111(int i,int j,int k,int ii,int jj,int kk)
{
	double kT = 0.000086173*temperature;
	double delta = em_bias;		// NEEDS TO BE FIGURED OUT
								// global variable em_bias is input value
	double E = 0.0;
	
	// Hard code EM bias for nearest neighbors
	
	// In-plane neighbors (bias is same for all layers)
	if (kk == k) {
		if (jj == j) {
			if (ii == i+1) E = delta;					// (i+1,j,k)
			else if (ii == i-1) E = (-1.0)*delta;		// (i-1,j,k)
		}
		else {
			if (ii == i+1) E = delta/2.0;				// (i+1,j-1,k),(i+1,j+1,k)
			else if (ii == i) E = (-1.0)*delta/2.0;		// (i,j-1,k),(i,j+1,k)
		}
	}
	else {	// Out-of-plane neighbors (bias depends on which layer)
		if (k % 3 == 0)			// A layer
		{
			if ( (ii == i   && jj == j   && kk == k+1) ||
				 (ii == i+1 && jj == j-1 && kk == k-1) ) E = delta/2.0;
			if ( (ii == i-1 && jj == j   && kk == k+1) ||
				 (ii == i   && jj == j-1 && kk == k-1) ) E = (-1.0)*delta/2.0;
		}
		else if (k % 3 == 1) 	// B layer
		{
			if (jj == j) {
				if (ii == i+1) E = delta/2.0;
				if (ii == i) E = (-1.0)*delta/2.0;
			}
		}
		else if (k % 3 == 2)	// C layer
		{
			if ( (ii == i   && jj == j   && kk == k-1) ||
				 (ii == i+1 && jj == j+1 && kk == k+1) ) E = delta/2.0;
			if ( (ii == i-1 && jj == j   && kk == k-1) ||
				 (ii == i   && jj == j+1 && kk == k+1) ) E = (-1.0)*delta/2.0;
		}
	}

	return exp(E/kT);
}

/* ----------------------------------------------------------------------
   Convert site (i,j,k) to physical coordinates (x,y,z) on fcc(111)
------------------------------------------------------------------------- */
void AppMigration::get_xyz_fcc_12nn_111(int i, int j, int k, 
										double *x, double *y, double *z)
{
	// global variable latt_const is the material's lattice constant, a

	double sq = 1.0 / sqrt(2.0);
	double sq3 = 1.0 / sqrt(3.0);

	// global array positions
	int gi = i + nx_offset;
	int gj = j + ny_offset;
	int gk = k + nz_offset;

	// x and y values depend on the gk layer
	if (gk % 3 == 0)		// A layer
	{
		// x value also depends on even or odd row (gj)
		if (gj % 2 == 0) *x = (double)(gi)*sq*latt_const;
		else *x = ((double)(gi)+0.5)*sq*latt_const;

		*y = (double)(gj)*latt_const/2.0;
	}
	else if (gk % 3 == 1)	// B layer
	{
		// x value also depends on even or odd row (gj)
		if (gj % 2 == 0) *x = ((double)(gi)+0.5)*sq*latt_const;
		else *x = ((double)(gi)+1.0)*sq*latt_const;

		*y = ((double)(gj)/2.0 + (1.0/6.0))*latt_const;
	}
	else if (gk % 3 == 2) 	// C layer
	{
		// x value also depends on even or odd row (gj)
		if (gj % 2 == 0) *x = (double)(gi)*sq*latt_const;
		else *x = ((double)(gi)+0.5)*sq*latt_const;

		*y = ((double)(gj)/2.0 + (1.0/3.0))*latt_const;
	}

	// z value is the same for all gi, gj
	*z = (double)(gk)*latt_const*sq3;
}

/* ----------------------------------------------------------------------
   Calculate box upper bounds for visualization on fcc(111)
------------------------------------------------------------------------- */
void AppMigration::box_bounds_fcc_12nn_111(double *xhi,double *yhi,double *zhi)
{
	// global variable latt_const is the material's lattice constant, a

	double sq = 1.0 / sqrt(2.0);
	double sq3 = 1.0 / sqrt(3.0);

	*xhi = ((double)(nx_global) + 0.5)*sq*latt_const;
	*yhi = (((double)(ny_global)*0.5)+(1.0/3.0))*latt_const;
	*zhi = (double)(nz_global)*sq3*latt_const;
}

/* ----------------------------------------------------------------------
   Get nearest neighbor sites for site (i,j,k) on hcp(0001)
   (0001) surface normal to z axis, (1-100) along x axis
------------------------------------------------------------------------- */
void AppMigration::find_nnsites_hcp_12nn_0001(int *nnsites, int i, int j, int k)
{
	int n, ni[12], nj[12], nk[12];

	// Hard code indices for nearest neighbors
	if (k % 2 == 0)		// even layer
	{
		ni[0] = i-1; nj[0] = j-1; nk[0] = k-1;
		ni[1] = i-1; nj[1] = j-1; nk[1] = k  ;
		ni[2] = i-1; nj[2] = j-1; nk[2] = k+1;
		ni[3] = i-1; nj[3] = j  ; nk[3] = k-1;
		ni[4] = i-1; nj[4] = j  ; nk[4] = k  ;
		ni[5] = i-1; nj[5] = j  ; nk[5] = k+1;
		ni[6] = i  ; nj[6] = j-1; nk[6] = k  ;
		ni[7] = i  ; nj[7] = j  ; nk[7] = k-1;
		ni[8] = i  ; nj[8] = j  ; nk[8] = k+1;
		ni[9] = i  ; nj[9] = j+1; nk[9] = k  ;
		ni[10]= i+1; nj[10]= j-1; nk[10]= k  ;
		ni[11]= i+1; nj[11]= j  ; nk[11]= k  ;
	}
	else				// odd layer
	{
		ni[0] = i-1; nj[0] = j-1; nk[0] = k  ;
		ni[1] = i-1; nj[1] = j  ; nk[1] = k  ;
		ni[2] = i  ; nj[2] = j-1; nk[2] = k  ;
		ni[3] = i  ; nj[3] = j  ; nk[3] = k-1;
		ni[4] = i  ; nj[4] = j  ; nk[4] = k+1;
		ni[5] = i  ; nj[5] = j+1; nk[5] = k  ;
		ni[6] = i+1; nj[6] = j-1; nk[6] = k-1;
		ni[7] = i+1; nj[7] = j-1; nk[7] = k  ;
		ni[8] = i+1; nj[8] = j-1; nk[8] = k+1;
		ni[9] = i+1; nj[9] = j  ; nk[9] = k-1;
		ni[10]= i+1; nj[10]= j  ; nk[10]= k  ;
		ni[11]= i+1; nj[11]= j  ; nk[11]= k+1;
	}

	// Convert neighbor indices to sites
	for (n = 0; n < 12; n++) {
		nnsites[n] = ijk2gsite(ni[n],nj[n],nk[n]);
	}
}

/* ----------------------------------------------------------------------
   Calculate EM bias on nearest neighbor sites of (i,j,k) on hcp(0001)
------------------------------------------------------------------------- */
double AppMigration::calc_bias_hcp_12nn_0001(int i,int j,int k,int ii,int jj,int kk)
{
	double kT = 0.000086173*temperature;
	double delta = em_bias;		// NEEDS TO BE FIGURED OUT
								// global variable em_bias is input value
	double E = 0.0;
	
	// Hard code EM bias for nearest neighbors
	
	// In-plane neighbors (bias is same for all layers)
	if (kk == k) {
		if (ii == i+1) E = delta;
		else if (ii == i-1) E = (-1.0)*delta;
	}
	else {	// Out-of-plane neighbors (bias depends on which layer)
		if (k % 2 == 0)			// even layer
		{
			if (ii == i) E = 0.5*delta;
			else if (ii == i-1) E = (-0.5)*delta;
		}
		else					// odd layer
		{
			if (ii == i+1) E = 0.5*delta;
			else if (ii == i) E = (-0.5)*delta;
		}
	}

	return exp(E/kT);
}

/* ----------------------------------------------------------------------
   Convert site (i,j,k) to physical coordinates (x,y,z) on hcp(0001)
------------------------------------------------------------------------- */
void AppMigration::get_xyz_hcp_12nn_0001(int i, int j, int k, 
										 double *x, double *y, double *z)
{
	// global variable latt_const is the material's lattice constant, a
	// global variable c_to_a is the ratio c/a for hexagonal lattices
	// c_to_a HAS TO BE ADDED IN CONSTRUCTOR, HEADER FILE
	// c = lattice spacing in z direction, a = lattice spacing in x and y 

	double sq = 1.0 / sqrt(2.0);
	double c_to_a = 1.0;
	double c = c_to_a*latt_const;

	// global array indices
	int gi = i + nx_offset;
	int gj = j + ny_offset;
	int gk = k + nz_offset;
	
	// x value depends on even or odd layer (gk)
	if (gk % 2 == 0) *x = (double)(gi)*sq*latt_const;	// even layer
	else *x = ((double)(gi)+0.5)*sq*latt_const;			// odd layer

	// y value depends on even or odd column (gi)
	if (gi % 2 == 0) *y = ((double)(gj)+0.5)*latt_const;// even column
	else *y = (double)(gj)*latt_const;					// odd column

	// z value is the same for all gi, gj
	*z = (double)(gk)*0.5*c;	
}

/* ----------------------------------------------------------------------
   Calculate box upper bounds for visualization on hcp(0001)
------------------------------------------------------------------------- */
void AppMigration::box_bounds_hcp_12nn_0001(double *xhi,double *yhi,double *zhi)
{
	// global variable latt_const is the material's lattice constant, a
	// global variable c_to_a is the ratio c/a for hexagonal lattices
	// c_to_a HAS TO BE ADDED IN CONSTRUCTOR, HEADER FILE
	// c = lattice spacing in z direction, a = lattice spacing in x and y 

	double sq = 1.0 / sqrt(2.0);
	double c_to_a = 1.0;
	double c = c_to_a*latt_const;

	*xhi = ((double)(nx_global) + 0.5)*sq*latt_const;
	*yhi = ((double)(ny_global) + 0.5)*latt_const;
	*zhi = (double)(nz_global)*0.5*c;
}

/* ---------------------------------------------------------------------- */

void AppMigration::xyz(int i, int j, int k, double *x, double *y, double *z)
{
  (this->*xyz_fp)(i,j,k,x,y,z);
}

/* ---------------------------------------------------------------------- */

void AppMigration::box_bounds(double *xlo, double *xhi, double *ylo,
			      double *yhi, double *zlo, double *zhi)
{
  *xlo = 0.0;
  *ylo = 0.0;
  *zlo = 0.0;

  (this->*box_fp)(xhi,yhi,zhi);
}

/* ---------------------------------------------------------------------- */
