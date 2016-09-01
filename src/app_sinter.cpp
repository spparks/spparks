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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_sinter.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "output.h"

#include <cfloat>
#include <cassert>
#include <cctype>
#include <new>
#include <fstream>

#include <iostream>
using std::cout;

using namespace SPPARKS_NS;

#define ONE_POSITIVE(a,b) ( a > 0 ? a : b )

// Taken from Class Lattice file lattice.cpp
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};
// Taken from Class AppLattice file app_lattice.cpp	   
enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};	   


/* ---------------------------------------------------------------------- */

AppSinter::AppSinter(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 0;
  delpropensity = 2; // Leave delpropensity and delevent like that to avoid producing the "Ghost connection was not found" error
  delevent = 1;
  allow_kmc = 0;
  allow_rejection = 1;
  allow_masking = 0; // ??
  numrandom = 1; // ??

  create_arrays();

  double gg_events=1.0, pm_events=1.0, a_events=1.0;
  double total_events = gg_events + pm_events + a_events;	  	

  gg_frequency = gg_events / total_events;  
  pm_frequency = pm_events / total_events;
  a_frequency = a_events / total_events;
  gb_factor = 1.0;
  gb_ini = -1.0; // Force to calculate it when starting sintering 
  dt_sweep = 1.0 / gg_frequency;

  vac_made = 0;
  
  frame_depth = 1;
  border_depth = 2;

  temperature = pore_migration_temperature = 1.0;
  annihilation_temperature = 15.0 * pore_migration_temperature;  // this is somewhat arbitrary, anni_temp can be independent of PM_temp.***********

  t_inverse = 1.0 / temperature;
  t_pm_inverse = 1.0 / pore_migration_temperature;
  t_anni_inverse = 1.0 / annihilation_temperature;
  time_sinter_start = 50;
  
  size_annihilist = 0;
  size_collapsinglist = 0;
}

/* ---------------------------------------------------------------------- */

AppSinter::~AppSinter()
{

  unique.clear();
  
  adjacent_spins.clear();
  grain_start.clear();
  if ( nprocs > 1 ) {
		check.clear();
		annihilation_list.clear();
		annihilation_spin.clear();
		annihilist_dist.clear();
		displacement.clear();
	}		
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppSinter::grow_app()
{
  spin = iarray[0];
}

/* ---------------------------------------------------------------------- */

void AppSinter::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"event_ratios") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal event_ratios command (provide three positive numbers for gg, pm and annihilation events respectively)");
	 double gg_events = atof( arg[0] );  // If not in inputfile, then default values provided before.
	 double pm_events = atof( arg[1] );	
	 double a_events = atof( arg[2] );
	 
	 if ( gg_events < 0 ) error->all(FLERR,"Illegal event_ratios command (number of gg events cannot be negative)");
	 if ( pm_events < 0 ) error->all(FLERR,"Illegal event_ratios command (number of pm events cannot be negative)");
 	 if ( a_events < 0 ) error->all(FLERR,"Illegal event_ratios command (number of annihilation events cannot be negative)");
		 	 
	 double total_events = gg_events + pm_events + a_events;	  	
	 gg_frequency = gg_events / total_events;  
	 pm_frequency = pm_events / total_events;
	 a_frequency = a_events / total_events;
	 gb_factor = 1.0;	 
	 dt_sweep = 1.0 / gg_frequency;
  }
  else if (strcmp(command,"events_temperatures") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal events_temperatures command (provide three positive numbers for gg, pm and annihilation events respectively)");
	 temperature = atof( arg[0] );
	 pore_migration_temperature = atof( arg[1] );
	 annihilation_temperature = atof( arg[2] );	 	

	 if ( temperature < 0 ) error->all(FLERR,"Illegal events_temperatures command (temperature for gg cannot be negative)");
	 if ( pore_migration_temperature < 0 ) error->all(FLERR,"Illegal events_temperatures command (temperature for pm cannot be negative)");
	 if ( annihilation_temperature < 0 ) error->all(FLERR,"Illegal events_temperatures command (temperature for annihilation cannot be negative)");
	 
	 t_inverse = 1.0 / temperature;
    t_pm_inverse = 1.0 / pore_migration_temperature;
    t_anni_inverse = 1.0 / annihilation_temperature;
  }
   else if (strcmp(command,"time_sinter_start") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal time_sinter_start command (provide one positive number)");
	 time_sinter_start = atof( arg[0] );	 
	 if ( time_sinter_start < 0 ) error->all(FLERR,"Illegal time_sinter_start command (cannot be negative)");
  } else error->all(FLERR,"Unrecognized command");
  
}

/* ---------------------------------------------------------------------- */

void AppSinter::init_app()
{

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;
	try {
		unique.resize( 1 + maxneigh ); // making sure that the vectors have sufficient memory.
		adjacent_spins.resize( 1 + maxneigh );
		grain_start.resize( 1 + maxneigh );
		if ( nprocs > 1 ) {
			check.resize( nlocal+nghost );
			
      	// MAX is a macro defined in 'pointers.h'
			int max_side = MAX(nx, MAX(ny, nz) );   
			const int SIZE_S ( max_side * 100 );

			annihilation_list.resize( SIZE_S );
			annihilation_spin.resize( SIZE_S );
			annihilist_dist.resize ( nprocs );
			displacement.resize( nprocs );
//			collapsing_directions.capacity( SIZE_S*1000 );
			collapsing_directions.resize( SIZE_S*8 );
			
			for (int i = 0; i < nlocal; i++)
				hash.insert(std::pair<int,int> (id[i],i));

		}
		dt_sweep = 1.0 / gg_frequency;
//		overimpose_frame();
	}	
	catch (std::bad_alloc& ba) {
		error->one(FLERR,"bad_alloc caught. Failed to allocate array in app_sinter init_app");
	}

}

/* ---------------------------------------------------------------------- */

void AppSinter::setup_app()  // making sure that some i may be listed for annihilation
{

  	dimension = domain->dimension;
	nbasis = domain->lattice->nbasis;
  
   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

  	Dx = (domain->boxxhi - domain->boxxlo) / nx; // Cristina thought this may become important in the future *****
	Dy = (domain->boxyhi - domain->boxylo) / ny; // so unless designated in the inputfile, Dx = 1/nx, Dy = 1/ny ....
	Dz = (domain->boxzhi - domain->boxzlo) / nz;
	
	assert(Dx > 0);  // for error checking
	if ( dimension > 1 )	assert(Dy > 0);
	if ( dimension > 2 )	assert(Dz > 0);
 
	nx_procs = domain->procgrid[0];
	ny_procs = domain->procgrid[1];
	nz_procs = domain->procgrid[2];
 
  	xgrid_proc = (double) nx / nx_procs; // Cristina, what is nx, ny, nz and where it is defined? nx = size of sim in X-dir.
	ygrid_proc = (double) ny / ny_procs; // xgrid_proc is the size in X at each processor.
	zgrid_proc = (double) nz / nz_procs; // nx_procs defined in create_lattice (), but not xgrid_proc.
	
	if ( nprocs > 1 ) {
		for (int i = 0; i < nlocal+nghost; i++) {
			check[i] = false;
		}	
//if( me == 0 ) printf("nx: %d ny: %d nz: %d, Dx: %lf Dy: %lf Dz: %lf\n", nx, ny, nz, Dx, Dy, Dz);
//printf("in proc: %d my position in proc grid is: %d %d %d\n", me, domain->myloc[0], domain->myloc[1], domain->myloc[2]);		
//printf("in proc: %d my portion of box is: x: %lf %lf, y: %lf %lf, z: %lf %lf\n", me, domain->subxlo, domain->subxhi, domain->subylo, domain->subyhi, domain->subzlo, domain->subzhi);
//printf("check list has been initialized to false\n");
	}	
	
	gb_ini = -1.0; // Force to calculate it when starting sintering 
	initialize_parameters_density_calculation();
//	check_space_initialization();
	overimpose_frame();
}

/* ---------------------------------------------------------------------- */

void AppSinter::initialize_parameters_density_calculation()
{
	// Determine the central parallelepiped for calculating density
/*	double gs = count_grain_sites();
	double occupied_fraction = (double) gs / (double)( nx * ny * nz );
	double rcube_fraction = pow( occupied_fraction, 1. / 3.);
	int nx_density = (int)floor( nx * rcube_fraction );
	int ny_density = (int)floor( ny * rcube_fraction );
	int nz_density = (int)floor( nz * rcube_fraction );
  
	// Open interval: x > xstart_density and x < xend_density ...
	xstart_density = (nx - nx_density) / 2;
	xend_density = xstart_density + nx_density;
	ystart_density = (ny - ny_density) / 2;
	yend_density = ystart_density + ny_density;
	zstart_density = (nz - nz_density) / 2;
	zend_density = zstart_density + nz_density; */

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

  int xsize = nx / 3;
  int ysize = ny / 3;
  int zsize = nz / 3;
  
  xstart_density = nx*0.33;
  xend_density = xstart_density + xsize;
  ystart_density = ny*0.33;
  yend_density = ystart_density + ysize;
  zstart_density = nz*0.33;
  zend_density = zstart_density + zsize; 

}

/* ---------------------------------------------------------------------- */

void AppSinter::overimpose_frame()
{	
	int total_border_depth = frame_depth + border_depth;
	int latstyle = domain->lattice->style;

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	int i,j,k;
	int isite;		
	int flag = 0;
	for (int ind = 0; ind < nlocal; ind++) {
		global_to_grid( id[ind], i, j, k );
		if (latstyle == LINE_2N) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ) {
				spin[ind] = FRAME;
			}	
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ) {
				spin[ind] = VACANT;
			}	
		} 
		else if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ||
				(j < frame_depth) || (abs(j-ny+1) < frame_depth) ) {
				spin[ind] = FRAME;
			}	
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ||
					(j < total_border_depth) || (abs(j-ny+1) < total_border_depth) ) {
				spin[ind] = VACANT;
			}
		}	
		else if (latstyle == SC_6N || latstyle == SC_26N || 
				latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ||
				(j < frame_depth) || (abs(j-ny+1) < frame_depth) ||
				(k < frame_depth) || (abs(k-nz+1) < frame_depth) )	{
				spin[ind] = FRAME;
			}
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ||
					(j < total_border_depth) || (abs(j-ny+1) < total_border_depth) ||
					(k < total_border_depth) || (abs(k-nz+1) < total_border_depth) ) {
				spin[ind] = VACANT;
			}
		}		
	}
}

/* ----------------------------------------------------------------------
	Check space for frame, border and actual simulation 
	All the unit cell is considered as part of the frame / border, thus
	for lattice structures with more than one base node all the nodes in
	the same unit cell are initialized equally at the frame / border and
	randomly (different) for any other location
------------------------------------------------------------------------- */

void AppSinter::check_space_initialization()
{	
	int total_border_depth = frame_depth + border_depth;
	int latstyle = domain->lattice->style;

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

//printf("Inside check_space_initialization\n");
//printf("nglobal: %d nlocal: %d\n", nglobal, nlocal );

	int i,j,k;
	int isite;		
	int flag = 0;
	for (int ind = 0; ind < nlocal; ind++) {
//printf("%d ", i);	
		isite = spin[ind];
		global_to_grid( id[ind], i, j, k );
//printf("i:%d j:%d k:%d\n", i, j, k );
		if (latstyle == LINE_2N) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ) {
				if ( isite != FRAME ) {
					flag = 1;
					break;
				}	
			}	
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ) {
				if ( isite != VACANT ) {
					flag = 1;
					break;
				}	
			}	
			else {
				if (isite < VACANT) {
					flag = 1;
					break;
				}	
			}	
		} 
		else if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ||
				(j < frame_depth) || (abs(j-ny+1) < frame_depth) ) {
				if ( isite != FRAME ) {
					flag = 1;
					break;
				}
			}	
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ||
					(j < total_border_depth) || (abs(j-ny+1) < total_border_depth) ) {
				if ( isite != VACANT ) {
					flag = 1;
					break;
				}	
			}
			else {
				if (isite < VACANT) {
					flag = 1;
					break;
				}	
			}	
		}	
		else if (latstyle == SC_6N || latstyle == SC_26N || 
				latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
			if ( (i < frame_depth) || (abs(i-nx+1) < frame_depth) ||
				(j < frame_depth) || (abs(j-ny+1) < frame_depth) ||
				(k < frame_depth) || (abs(k-nz+1) < frame_depth) )	{
				if ( isite != FRAME ) {
					flag = 1;
					break;
				}	
			}
			else if ( (i < total_border_depth) || (abs(i-nx+1) < total_border_depth) ||
					(j < total_border_depth) || (abs(j-ny+1) < total_border_depth) ||
					(k < total_border_depth) || (abs(k-nz+1) < total_border_depth) ) {
				if ( isite != VACANT ) {
					flag = 1;
					break;
				}	
			}
						
			else {
				if (isite < VACANT) {
					flag = 1;
					break;
				}	
			}	
		}		
	}
	
	int flagall;
	MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
	if (flagall) error->all(FLERR,"One or more sites have invalid values");
}


/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppSinter::site_energy(int i)
{
  int isite = spin[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != spin[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
   propensity for one event is based on einitial,efinal
------------------------------------------------------------------------- */
double AppSinter::site_propensity(int i)
{
  return 0;
}


/* ----------------------------------------------------------------------
   choose and perform an exchange for site
   update propensities of all affected sites
   ignore neighbor sites that should not be updated (isite < 0)
------------------------------------------------------------------------- */

void AppSinter::site_event(int i, class RandomPark *random)
{
}

/* ----------------------------------------------------------------------
   perform a site event with rejection
------------------------------------------------------------------------- */

void AppSinter::site_event_rejection(int i, RandomPark *random)
{

	// event = exchange with random neighbor
	int oldstate = spin[i];
	
	if (oldstate == FRAME) 		return;
	
	if ( nprocs > 1 && check[i] ) return; // Vacancy being annihilated
	
	if(oldstate == VACANT) { // ****ASK, is i the vacancy/pore site always? **** 
		
		if ( (time > time_sinter_start) && (random->uniform() < a_frequency*gb_factor) ) {
//printf("Proc %d is starting to make vacany in position %d (local) %d (global)\n", me, i, id[i]);			
			make_vacancy(i, random);
		}		
		else if (random->uniform() < pm_frequency) { // Pore Migration step
			// Select random neighbor
			int iran = (int) (numneigh[i]*random->uniform());
			if (iran >= numneigh[i]) iran = numneigh[i] - 1;
			int j = neighbor[i][iran];
			int neighstate = spin[j];
		
			if(neighstate > VACANT) { // If the neighbor selected is a grain site
				double einitial = site_energy(i) + site_energy(j);
				spin[i] = neighstate;
				spin[j] = oldstate;
				int gs = i; // keep track of the new position for grain site
				
				// Initialize minimum energy with value of energy for 
				// new gran site (after making the exchange)
				double min_energy = site_energy(gs); 
				int newstate;
				choose_neighbor_grain_site_minimizing_energy( gs, neighstate, random, min_energy, newstate );
				spin[gs] = newstate;
				double efinal = min_energy + site_energy(j); //site_energy(i) + site_energy(j);

				// accept or reject the event
				if (efinal <= einitial) {
				} 
				else if (pore_migration_temperature == 0.0) {
					spin[i] = oldstate;
					spin[j] = neighstate;
				} 
				else if (random->uniform() > exp((einitial-efinal)*t_pm_inverse)) {
					spin[i] = oldstate;
					spin[j] = neighstate;
				}
			}
		}
	}
	else if (random->uniform() < gg_frequency) { // if oldstate > 1 i.e. if it corresponds to a grain site	
		// Calculate a GG step
		double einitial = site_energy(i); 
		
		int nevent = 0, m;
		for(int j = 0; j < numneigh[i]; j++) {
			int value = spin[neighbor[i][j]];
			if ( value == FRAME ) continue; // Outside simulation space
			if ( value == VACANT || value == oldstate ) continue;		// Vacancy or Same spin
			for ( m = 0; m < nevent; m++ )
				if ( value == unique[m] ) break;		// Registered event
			if ( m < nevent ) continue;	// Previous cycle was interrupted because event was found
			unique[nevent++] = value;
		}
      int iran = (int) (nevent*random->uniform());
      if(iran >= nevent) return;//iran = nevent - 1; --> protection for nevent = 0
      spin[i] = unique[iran];
      double efinal = site_energy(i);

      if(efinal <= einitial){
      }
      else if (temperature == 0.0) 
			spin[i] = oldstate;
      else if (random->uniform() > exp((einitial-efinal)*t_inverse)) 
			spin[i] = oldstate;
	}
	
	if ( spin[i] != oldstate ) naccept++;
}
/* ----------------------------------------------------------------------
 driver of rejection kinetic Monte Carlo
 ------------------------------------------------------------------------- */ 
void AppSinter::iterate_rejection(double stoptime)
{ 
/*  app_lattice has functions/methods that app_sinter uses.  However, it need to stop simulation
   *  at the end of each sector, which app_lattice does not.  So, Cristina told app_lattice to expect
   *  this to be added/modified in app_sinter.  And here it is.  Simulations at each processor stops,
   *  performs annihilations and restarts.
   *
   *  This function actually performs everything - annihilation, pm, gg
   */
   
	int i,icolor,nselect,nrange,jset;
	int *site2i;
	
	// set loop is over:
	// sectors if there are sectors and no colors
	// colors if there are colors and no sectors
	// first nsector sets if there are both sectors and colors
	
	int nset_loop = nset;
	if (bothflag) nset_loop = nsector;
	
	double nextgbupdate = time_sinter_start;
	
	int done = 0;
	while (!done) {
		for (int iset = 0; iset < nset_loop; iset++) {
			if (nprocs > 1) {
				timer->stamp();
				if (sectorflag) comm->sector(iset);
				else comm->all();
				timer->stamp(TIME_COMM);
			}
			
			if (Lmask) boundary_clear_mask(iset);
			
			timer->stamp();
			
			// sectors but no colors (could also be no sectors)
			// random selection of sites in iset
			
			if (sweepflag == RANDOM) {
				site2i = set[iset].site2i;
				nrange = set[iset].nlocal;
				nselect = set[iset].nselect;
				for (i = 0; i < nselect; i++) 
					sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
				(this->*sweep)(nselect,sitelist);
				nattempt += nselect;

				// sectors but no colors, or colors but no sectors
				// ordered sweep over all sites in iset

			} else if (bothflag == 0) {
				for (i = 0; i < set[iset].nloop; i++)
					(this->*sweep)(set[iset].nlocal,set[iset].site2i);
				nattempt += set[iset].nselect;

				// sectors and colors
				// icolor loop is over all colors in a sector
				// jset = set that contains sites of one color in one sector
				// ordered sweep over all sites in jset

			} else {
				for (icolor = 0; icolor < ncolors; icolor++) {
					jset = nsector + iset*ncolors + icolor;
					for (i = 0; i < set[jset].nloop; i++)
						(this->*sweep)(set[jset].nlocal,set[jset].site2i);
					nattempt += set[jset].nselect;
			}
		}
		
		timer->stamp(TIME_SOLVE);
		if (nprocs > 1) {
			if (sectorflag) comm->reverse_sector(iset);
			else comm->all_reverse();
			timer->stamp(TIME_COMM);
			update_annihilations(ranapp);
			timer->stamp(TIME_APP);
		}
	}

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
	if ( time > nextgbupdate ) nextgbupdate = calculate_gb_update( time );
  }   
}


/////////////////////////////////////////// MAPPINGS ///////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Mappings between real space, simulation grid and processor grid
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
	Given a set of coordinates this function returns the processor that
	owns the closest lattice site (always mapped to the base point of the 
	unit cell (0,0,0) ) and the integer coordinates of the site in the
	lattice 
------------------------------------------------------------------------- */

int AppSinter::xyz_to_processor_and_grid( double x, double y, double z, int &grid_x, int &grid_y, int &grid_z )
{

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	grid_x = (int)floor((x-domain->boxxlo) / Dx + 0.5);
	grid_y = (int)floor((y-domain->boxylo) / Dy + 0.5);
	grid_z = (int)floor((z-domain->boxzlo) / Dz + 0.5);
//	if ( grid_x < 0 ) grid_x = 0;
//	if ( grid_y < 0 ) grid_y = 0;
//	if ( grid_z < 0 ) grid_z = 0;
	
	if ( grid_x >= nx ) grid_x = nx-1;
	if ( grid_y >= ny ) grid_y = ny-1;
	if ( grid_z >= nz ) grid_z = nz-1;
	
	int procx = (int)(grid_x / xgrid_proc);
	if ( procx > nx_procs ) procx = nx_procs-1;
	int procy = (int)(grid_y / ygrid_proc);
	if ( procy > ny_procs ) procy = ny_procs-1;
	int procz = (int)(grid_z / zgrid_proc);
	if ( procz > nz_procs ) procz = nz_procs-1;
		
	//int proc = procx * ny_procs * nz_procs + procy * nz_procs + procz; // OLD implementation of SPPARKS
	int proc = procz * nx_procs * ny_procs + procy * nx_procs + procx;
	
	assert ( proc >= 0 && proc < nprocs);
	
	return proc;

} 

/* ----------------------------------------------------------------------
	Given a set of coordinates this function returns the processor that
	owns the closest lattice site (always mapped to the base point of the 
	unit cell (0,0,0) ) and the global index of the site
------------------------------------------------------------------------- */

int AppSinter::xyz_to_processor_and_global( double x, double y, double z, int &iglobal )
{

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	int grid_x = (int)floor((x-domain->boxxlo) / Dx + 0.5);
	int grid_y = (int)floor((y-domain->boxylo) / Dy + 0.5);
	int grid_z = (int)floor((z-domain->boxzlo) / Dz + 0.5);
//	if ( grid_x < 0 ) grid_x = 0;
//	if ( grid_y < 0 ) grid_y = 0;
//	if ( grid_z < 0 ) grid_z = 0;
	
	if ( grid_x >= nx ) grid_x = nx-1;
	if ( grid_y >= ny ) grid_y = ny-1;
	if ( grid_z >= nz ) grid_z = nz-1;
	
	int procx = (int)(grid_x / xgrid_proc);
	if ( procx > nx_procs ) procx = nx_procs-1;
	int procy = (int)(grid_y / ygrid_proc);
	if ( procy > ny_procs ) procy = ny_procs-1;
	int procz = (int)(grid_z / zgrid_proc);
	if ( procz > nz_procs ) procz = nz_procs-1;
	
//	int proc = procx * ny_procs * nz_procs + procy * nz_procs + procz; // OLD implementation of SPPARKS
	int proc = procz * nx_procs * ny_procs + procy * nx_procs + procx;

	assert ( proc >= 0 && proc < nprocs);
	
	int ind_m = 0; // To take into acount other types of lattices but always mapping to base (0,0,0)
	int int_new_index = grid_z * nx * ny * nbasis + grid_y * nx * nbasis + grid_x * nbasis + ind_m; //nbasis in SPPARKS is nbasis of lattice theory.
	
	assert ( int_new_index >= 0 && int_new_index < nglobal);
	iglobal = int_new_index + 1;
	
//printf ("In proc: %d I calculated that iglobal: %d is in proc: %d, from owner array it is proc: %d\n", me, iglobal, proc, owner[iglobal]);
	
	return proc;
} 

/* ----------------------------------------------------------------------
   Given a set of coordinates this function returns the local index of 
	the closest lattice site (always mapped to the base point of the 
	unit cell (0,0,0) )
	It fails if the point is not inside the local box
------------------------------------------------------------------------- */

int AppSinter::xyz_to_local( double x, double y, double z )
{
   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	int grid_x = (int)floor((x-domain->boxxlo) / Dx + 0.5);
	int grid_y = (int)floor((y-domain->boxylo) / Dy + 0.5);
	int grid_z = (int)floor((z-domain->boxzlo) / Dz + 0.5);
//	if ( grid_x < 0 ) grid_x = 0;
//	if ( grid_y < 0 ) grid_y = 0;
//	if ( grid_z < 0 ) grid_z = 0;
	
	if ( grid_x >= nx ) grid_x = nx-1;
	if ( grid_y >= ny ) grid_y = ny-1;
	if ( grid_z >= nz ) grid_z = nz-1;	
	
	
	int ind_m = 0; // To take into acount other types of lattices but always mapping to base (0,0,0)
		
	int int_new_index = grid_z * nx * ny * nbasis + grid_y * nx * nbasis + grid_x * nbasis + ind_m;
	
	assert ( int_new_index >= 0 && int_new_index < nglobal);

	if ( nprocs > 1 ) {

		int iglobal = int_new_index + 1;

		std::map<int,int>::iterator loc; // find the local index from global, 2nd place is loc index.
		loc = hash.find(iglobal);
	
		int local_ind = loc->second;
//printf("in proc %d calculated correspondance: iglobal-ilocal: %d-%d, registered: iglobal-ilocal: %d-%d\n", me, iglobal, local_ind, id[local_ind], local_ind);		

		assert( id[local_ind] == iglobal );
//		assert ( local_ind >= 0 && local_ind < nlocal);
//printf("in proc %d the input coordinates are: %f %f %f, I think this is close to site %d (local) %d (global) with coordinates: %f %f %f\n", me, x, y, z, local_ind, iglobal, xyz[local_ind][0], xyz[local_ind][1], xyz[local_ind][2]);		
	
		return local_ind;
	}
	return int_new_index;	
}

/* ----------------------------------------------------------------------
	Given the global index of the site this function returns the integer 
	coordinates of the site in the lattice 
------------------------------------------------------------------------- */

void AppSinter::global_to_grid( int iglobal, int & i, int & j, int & k )
{
   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

//cout << "Inside global to grid\n";
	int latstyle = domain->lattice->style;
	nbasis = domain->lattice->nbasis;

	if (latstyle == LINE_2N) {
		i = (iglobal-1)/nbasis % nx;
		j = 0; 
		k = 0;
	} else if (latstyle == SQ_4N || latstyle == SQ_8N || latstyle == TRI) {
		i = (iglobal-1)/nbasis % nx;
		j = (iglobal-1)/nbasis / nx;
		k = 0;
	}	else if (latstyle == SC_6N || latstyle == SC_26N || 
	     latstyle == FCC || latstyle == BCC || latstyle == DIAMOND) {
		i = (iglobal-1)/nbasis % nx;
		j = (iglobal-1)/nbasis / nx % ny;
		k = (iglobal-1)/nbasis / (nx*ny);
	}		
//cout << "The global id: " << iglobal << " corresponds to the following coordinates: " << i << ", " << j << ", " << k << "\n";	
}

////////////////////////////////////// GRAIN MINIMIZING ENERGY ////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Starting from a random position select the neighbor grain site spin
	that minimizes the energy of the site i
------------------------------------------------------------------------- */

void AppSinter::choose_neighbor_grain_site_minimizing_energy( int i, int neighstate, RandomPark *random, double & min_gs_energy, int & newstate )
{
	// Choose a random offset (every grain should have equal probability)
	int initnn = (int) (numneigh[i]*random->uniform());
	if (initnn >= numneigh[i]) initnn = numneigh[i] - 1;
	   
	// Go through all neighbors
	int nevent = 0, m;
	newstate = neighstate;
	for(int nbor=0; nbor < numneigh[i]; nbor++) {
		if ( min_gs_energy <= (numneigh[i] / 2) ) break; // No other neighbor can produce a lowest energy
		int gs_neigh = neighbor[i][(nbor+initnn)%numneigh[i]];
		int value = spin[gs_neigh];
		if ( value == FRAME ) continue; // Outside simulation space
		if ( value == VACANT || value == neighstate ) continue; // Vacancy or same spin
		for ( m = 0; m < nevent; m++ )
			if ( value == unique[m] ) break; // Registered Event
		if ( m < nevent ) continue; // Previous cycle was interrupted because event was found
		unique[nevent++] = value;
		spin[i] = value;
		double gs_eng = site_energy(i);
		if ( gs_eng < min_gs_energy ) {
			min_gs_energy = gs_eng;
			newstate = value;
		}	
	}
}

/////////////////////////////////////////// VACANCIES AND ANNIHILATIONS //////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   create a vacancy
------------------------------------------------------------------------- */

void AppSinter::make_vacancy(int i, RandomPark *random)
{
	int oldstate = spin[i];	// It is a pore site
	assert( oldstate == VACANT );
	
	// Select random neighbor
	int iran = (int) (numneigh[i]*random->uniform());
	if (iran >= numneigh[i]) iran = numneigh[i] - 1;
	int j = neighbor[i][iran];
	int neighstate = spin[j];
	
	if(neighstate > VACANT) { // If the neighbor selected is a grain site
		double einitial = site_energy(i) + site_energy(j);
		spin[i] = neighstate;
		spin[j] = oldstate;
		double new_gs_energy = site_energy(i);
		double new_pore_energy = site_energy(j);
		// Check that interchange produces a vacancy (isolated pore)
		int int_new_pore_energy = (int)floor( new_pore_energy );
		if ( int_new_pore_energy == numneigh[j] ) { // New site is a vacancy
			double min_gs_energy = new_gs_energy;
			int newstate;
			choose_neighbor_grain_site_minimizing_energy( i, neighstate, random, min_gs_energy, newstate );
			spin[i] = newstate;
			double efinal = min_gs_energy + new_pore_energy;

			// accept or reject the event

			if (efinal <= einitial) {
			} 
			else if (annihilation_temperature == 0.0) {
				spin[i] = oldstate;
				spin[j] = neighstate;
			} 
			else if (random->uniform() > exp((einitial-efinal)*t_anni_inverse)) {
				spin[i] = oldstate;
				spin[j] = neighstate;
			}
		}
		else { // New site is not a vacancy -> revert exchange
			spin[i] = oldstate;
			spin[j] = neighstate;
		}
	}
	
	if ( spin[i] != oldstate ) {
		naccept++;
		vac_made++;
		assert ( spin[j] == VACANT );
		annihilate(j, random);
	}	
}


/* ----------------------------------------------------------------------
   annihilate a vacancy
------------------------------------------------------------------------- */

void AppSinter::annihilate(int i, RandomPark *random)
{
//	int oldstate = lattice[i];	// It is a vacancy
	if ( nprocs > 1 ) {
		register_annihilation_event( i );
	}
	else {
		int adjacent_spin, adjacent_start;
		if ( vacancy_adjacent_grain( i, random, adjacent_spin, adjacent_start ) ) {
			// Calculate mass center for adjacent grain
			vector<double> coord_cm ( dimension );
			vector<double> new_pos( dimension );
			vector<double> p( dimension );
			double mint;
			
			double volume = calculate_mass_center_adjacent_grain(adjacent_start, adjacent_spin, coord_cm);
			// Calculate the new position for vacancy and shift all intermediate positions
			calculate_vacancy_new_position( i, coord_cm, p, mint, new_pos );		
			collapse_in_direction( i, p, mint, new_pos );
		}
		else {
			error->one(FLERR,"In annihilation, site to annihilate is a PORE not a VACANCY");
		}	
	}	
}

/* ----------------------------------------------------------------------
    In parallel execution annihilation cannot be completed until the
	end of the sweep, when all processors are waiting to communicate
	with each other. Thus, the index of the site with the vacancy that 
	is being annihilated is registered (i).
	Array check is used to pin the site i so no further changes can be
	done until the annihilation is processed. Index i can correspond to
	a local (owned) or a ghost site. 
------------------------------------------------------------------------- */

void AppSinter::register_annihilation_event(int i) 
{
	assert( spin[i] == VACANT ); // if false, then run is terminated with an error statement.
	assert( !check[i] );
	
	annihilation_list[size_annihilist++] = i;
	check[i] = true;
//printf("in proc %d the local size of annihilation list is %d and the new event added is for site: %d\n", me, size_annihilist, i);	
}

/* ----------------------------------------------------------------------
   In Parallel execution, consolidate the list of annihilations made by 
	each processor, calculate the center of mass for each adjacent grain,
	compute the new vacancy position and generate a list of collapsing
	paths. The collapsing list is calculated by the processor that owns
	the vacancy to be annihilated. The states in the path are filled
	by each processor after the list of sites has been broadcasted.
	If collapsing some path moves one of the vacancies registered to
	be annihilated, then this annihilation is removed from the list and
	no collapsing is computed.
------------------------------------------------------------------------- */
 
void AppSinter::update_annihilations(RandomPark *random)
{
	if ( time > time_sinter_start ) {
	
		try {	
			// check:	- The vaccancy to be annihilated is still a vacancy
			//				- The registered spin is still a surrounding spin if not -> re-calculate the adjacent grain
//printf("In proc %d the local size of the annihilation list is %d\n", me, size_annihilist);			
			generate_annihilation_spin_list(random);
		
			int size_annihilist_all=0;
			MPI_Allreduce(&size_annihilist,&size_annihilist_all,1,MPI_INT,MPI_SUM,world);
		
			if ( size_annihilist_all ) {	
/*							
				const int GROUP(dimension+1);
				vector<double> cm_buffer_out ( size_annihilist_all*GROUP );			
				
				calculate_mass_center_distributed_grains( size_annihilist_all, GROUP, cm_buffer_out ); //com associated with each vacancy is calculated. 
				generate_local_vacancy_collapsing_list( size_annihilist_all, GROUP, cm_buffer_out );
*/				
				generate_local_vacancy_collapsing_list( size_annihilist_all );
				filter_local_vacancy_collapsing_list(); // making sure that overlapping vacancies are filtered.
				
				int size_collapsinglist_all = 0;
				MPI_Allreduce(&size_collapsinglist,&size_collapsinglist_all,1,MPI_INT,MPI_SUM,world); // all processors know the total number of pending annihilations 
				if ( size_collapsinglist_all ) {
					collapse_pending_list ( size_collapsinglist_all );
				}	
				// Allow site to change after updating annihilation	
				for ( int i = 0; i < size_annihilist; i++ )	{
					check[annihilation_list[i]] = false;	
				}
				size_annihilist = 0;
			}
		}		
		catch (std::bad_alloc& ba) {
			error->all(FLERR,"bad_alloc caught inside update_annihiliations");
		}
	}			 
}

/////////////////////////////////////////// PROCESSING ANNIHILATIONS - SERIAL //////////////////////////////////////////////////

/* ----------------------------------------------------------------------
	Given the index i of a supposed vacancy the adjacent grain, defined
	as the grain with the spin that occupies the mayority of sites in the
	neighborhood of the vacancy, is determined.
	If the site is indeed a vacancy function returns true, if the site
	is not a vacancy because it is surrounded by other vacant neighbors
	then function returns false.
	The search for the adjacent grain starts from a random neighbor
	to avoid the bias in case of a tie. 
	The adjacent grain its specified through its spin and a position 
	occupied by it.	
------------------------------------------------------------------------- */

bool AppSinter::vacancy_adjacent_grain( int i, RandomPark *random, int & adjacent_grain_spin, int & adjacent_grain_start )
{	
	// Select random neighbor
	int iran = (int) (numneigh[i]*random->uniform());
	if (iran >= numneigh[i]) iran = numneigh[i] - 1;

	// Go through all neighbors starting from random neighbor
	int nevent = 0, m;
	for(int nbor=0; nbor < numneigh[i]; nbor++) {
		int gs_neigh = neighbor[i][(nbor+iran)%numneigh[i]];
		int value = spin[gs_neigh];
		if ( value == VACANT )	return false;	// Other vacant site in the neighboorhood -> it is not a vacancy
		if ( value == FRAME ) continue; // Outside simulation space
		for ( m = 0; m < nevent; m++ ) {
			if ( value == unique[m] ) {
				(adjacent_spins[m])++;
				break; // Registered Event
			}	
		}	
		if ( m < nevent ) continue; // Previous cycle was interrupted because event was found
		unique[nevent] = value;
		adjacent_spins[nevent] = 1;
		grain_start[nevent++] = gs_neigh;
	}
	if ( !nevent ) // No grain site found in the neighborhood -> it is not a vacancy
		return false;
		
	int adjacent_spin = unique[0];
	int num_neigh_adj = adjacent_spins[0];
	int adjacent_start = grain_start[0];
	for ( int adj = 1; adj < nevent; adj++ ) {
		if ( adjacent_spins[adj] > num_neigh_adj ) {
			adjacent_spin = unique[adj];
			num_neigh_adj = adjacent_spins[adj];
			adjacent_start = grain_start[adj];				
		}
	}
	adjacent_grain_spin = adjacent_spin;	
	adjacent_grain_start = adjacent_start;
	
	return true;
}


/* ----------------------------------------------------------------------
	For execution in ONE processor
	Calculate the mass center for the adjacent grain
	Index i corresponds to one site in adjacent grain
	A site is considered to be part of the adjacent grain if it has the
	same spin as the one given and if it belongs to the neighborhood of
	the initial site given, or to any of its neighbors
------------------------------------------------------------------------- */

double AppSinter::calculate_mass_center_adjacent_grain(int i, int isite, vector<double> & cm )
{
	
	vector<bool> registered ( nlocal, false ); 
	
	for ( int j = 0; j < dimension; j++ ) 
		cm[j] = 0.0;
		
	assert( isite != FRAME && isite != VACANT );	
		
	std::stack<int> grainstack;	
	grainstack.push(i);
	registered[i] = true;
	double vol = 0.0;
	
	int ii, k;
	
	while ( grainstack.size() ) {
		k = grainstack.top();
//		assert ( spin[k] == isite );
//		assert ( k >= 0 && k < nlocal );
		grainstack.pop();
		vol++;
		for ( int j = 0; j < dimension; j++ ) 
			cm[j] += xyz[k][j];
		for (int j = 0; j < numneigh[k]; j++) {
			ii = neighbor[k][j];
			if (spin[ii] == isite && registered[ii] == false) {
				grainstack.push(ii);
				registered[ii] = true;
			}
		}		
	}
	
	for ( int j = 0; j < dimension; j++ ) 
		cm[j] /= vol;	
		
	return vol;	
}

/* ----------------------------------------------------------------------
	Calculate the new position for the vacancy
	Index i corresponds to the current site of vacancy
	coord_cm corresponds to the center of mass of the adjacent grain
	Direction is always assumed to be from the vacancy current position
	to the center of mass of the adjacent grain to border of simulation
	space
	At the end: p is the direction vector, new_pos is vacancy new position
	and mint is the number of discrete steps to be taken from the old 
	to the new position 
------------------------------------------------------------------------- */

void AppSinter::calculate_vacancy_new_position( int i, vector<double> & cm, vector<double> & p, double &mint, vector<double> & new_pos )
{
	
	vector<double> old_pos ( dimension );
	for ( int dim = 0; dim < dimension; dim++ ) {
		old_pos[dim] = xyz[i][dim];
		p[dim] = cm[dim] - old_pos[dim];
	}
		
	double txhi, txlo, tyhi, tylo, tzhi, tzlo;	
	
	txhi = txlo = tyhi = tylo = tzhi = tzlo = DBL_MAX;	
	
	if ( fabs(p[0]) > DBL_MIN ) {
		txhi =  (domain->boxxhi - old_pos[0]) / p[0];
		txlo =  (domain->boxxlo - old_pos[0]) / p[0];
	}
	if ( fabs(p[1]) > DBL_MIN ) {	
		tyhi =  (domain->boxyhi - old_pos[1]) / p[1];
		tylo =  (domain->boxylo - old_pos[1]) / p[1];
	}
	if ( fabs(p[2]) > DBL_MIN ) {
		tzhi =  (domain->boxzhi - old_pos[2]) / p[2];
		tzlo =  (domain->boxzlo - old_pos[2]) / p[2];
	}	
	
	double mintx = ONE_POSITIVE(txhi,txlo);
	double minty = ONE_POSITIVE(tyhi,tylo);
	double mintz = ONE_POSITIVE(tzhi,tzlo);	
	
	// MIN is a macro defined in 'pointers.h'
	mint = MIN(mintx, MIN(minty, mintz));
	
	if ( mint > (0.7 *DBL_MAX) ) {
		printf("proc: %d vacancy: (%lf, %lf, %lf) cm: (%lf, %lf, %lf) p: (%lf, %lf, %lf)\n", me, old_pos[0], old_pos[1], old_pos[2], cm[0], cm[1], cm[2], p[0], p[1], p[2]);
		printf("proc: %d size_collapsinglist: %d , collapsing_directions.size(): %d\n", me, size_collapsinglist, (int)collapsing_directions.size() );
		error->one(FLERR,"Vacancy to be annihilated and center of mass too close. Increment time to start sintering");
	}	
	
	for ( int dim = 0; dim < dimension; dim++ ) 
		new_pos[dim] = old_pos[dim] + mint * p[dim];		
}

/* ----------------------------------------------------------------------
	For execution in ONE processor
	Move the vacancy to the new position shifting the intermediate 
	positions back from the vacancy new position to the old one.
	Given:	- Index of the current position of the vacancy (ind_vac)
				- Direction vector between old and new position (p)
				- Maximum multiplier for vector p to reach new position (limit_t)
				- Vacancy new position (new_pos) 
	Starting from the vacancy position calculate the number of discrete
	steps to reach the new position.
	Use direction vector p multiplied by the size of discrete step to 
	determine the next position in the direction of collapsing. Move
	the lattice state of that site to the previous site in the path.
	Stop when the border is found. Set the site before the border to
	VACANT value.
------------------------------------------------------------------------- */


void AppSinter::collapse_in_direction( int ind_vac, vector<double> & p, double limit_t, vector<double> & new_pos )
{
   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	int i_vac = ind_vac/nbasis % nx;
	int j_vac = ind_vac/nbasis / nx % ny;
	int k_vac = ind_vac/nbasis / (nx*ny);

	int ind_end = xyz_to_local( new_pos[0], new_pos[1], new_pos[2] );
	int i_end = ind_end/nbasis % nx;
	int j_end = ind_end/nbasis / nx % ny;
	int k_end = ind_end/nbasis / (nx*ny);
	
	assert( spin[ind_end] == FRAME );	
	
	p[0] = (xyz[ind_end][0] - xyz[ind_vac][0]) / limit_t;
	p[1] = (xyz[ind_end][1] - xyz[ind_vac][1]) / limit_t;
	p[2] = (xyz[ind_end][2] - xyz[ind_vac][2]) / limit_t;
	
	
	int stepsx = abs(i_end-i_vac);
	int stepsy = abs(j_end-j_vac);
	int stepsz = abs(k_end-k_vac);
	
	int steps = MAX( stepsx, MAX(stepsy, stepsz) );
	
	if ( steps < 2 ) // Just one step means to exchange vacancy and frame -> NO ! it destroys frame
		return;
	
	double delta_tt = limit_t / steps;
	
	double x, y, z;
	int current_index=ind_vac, next_index;
	x = xyz[ind_vac][0];
	y = xyz[ind_vac][1];
	z = xyz[ind_vac][2];
	
	for ( int i = 0; i < steps; i++ ) {
		x += delta_tt * p[0];
		y += delta_tt * p[1];
		z += delta_tt * p[2];
		next_index = xyz_to_local(x, y, z);
		assert( next_index != current_index );
		int before = spin[current_index];
		if ( spin[next_index] > FRAME && next_index != ind_end ) {
			spin[current_index] = spin[next_index];
		}
		else {
			spin[current_index] = VACANT;	
			break;
		}	
		current_index = next_index;	
	}
}

/////////////////////////////////////////// PROCESSING ANNIHILATIONS - PARALLEL /////////////////////////////////////////////////

/* ----------------------------------------------------------------------
	For PARALLEL execution
    In parallel execution annihilation cannot be completed until the
	end of the sweep, when all processors are waiting to communicate
	with each other. At the end of the sweep, all the vacancies
	being annihilated are checked to see if they are vacancies
	still. If so the adjacent grain is calculated
	and the corresponding spin is registered to be broadcasted later.
------------------------------------------------------------------------- */

void AppSinter::generate_annihilation_spin_list(RandomPark *random)
{
//printf("proc %d is in function generate_annihilation_spin_list\n", me );	
	vector<int> copy_annihilation_list( size_annihilist );
	for ( int i = 0; i < size_annihilist; i++ ) {
		copy_annihilation_list[i] = annihilation_list[i];
	} 
	
	int old_size_annihilist = size_annihilist;
	size_annihilist = 0;
	int new_neigh_spin, dum;
	// Disccard vacancies that are not longer there
	// Determine adjacent grain
	for ( int i = 0; i < old_size_annihilist; i++ ) {
		if ( vacancy_adjacent_grain( copy_annihilation_list[i], random, new_neigh_spin, dum ) ) {
			annihilation_list[size_annihilist] = copy_annihilation_list[i];
			annihilation_spin[size_annihilist++] = new_neigh_spin;
		}
		else { // It is no longer a vacancy
			naccept--;
			vac_made--;
			check[copy_annihilation_list[i]] = false;
		}	
	}
}


/* ----------------------------------------------------------------------
	For PARALLEL execution
	Given:	- Index of the current position of the vacancy (ind_vac)
				- Direction vector between old and new position (p)
				- Maximum multiplier for vector p to reach new position (limit_t)
				- Vacancy new position (new_pos) 
	
	If vacancy new position and old position are owner by the same
	processor: move vacancy to the new position shifting the intermediate 
	positions back from the vacancy new position to the old one. If 
	they are in different processors register in the collapsing list:
	initial position of vacancy, vector p, maximum multiplier, number
	of discrete steps, and owner of vacancy.
------------------------------------------------------------------------- */

void AppSinter::register_collapsing_event(int ind_vac, vector<double> & p, double limit_t, vector<double> & new_pos)
// Vacancy to be annihilated is in the processor entering this function 
{
	assert( spin[ind_vac] == VACANT );

   // num sites along each axis in lattice
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

	// From local index to global
	int global_ind_vac = id[ind_vac] - 1;
	
	// Indices of vacancy -> ind_vac needs to be a global index !
	int gridx_vac = global_ind_vac/nbasis % nx;
	int gridy_vac = global_ind_vac/nbasis / nx % ny;
	int gridz_vac = global_ind_vac/nbasis / (nx*ny);

	int gridx_end, gridy_end, gridz_end;
	int proc = xyz_to_processor_and_grid( new_pos[0], new_pos[1], new_pos[2], gridx_end, gridy_end, gridz_end );

	if ( proc == me && owner[ind_vac] == me ) {
//		int zeros=0;
//		for ( int i = 0; i < nlocal; ++i )
//			if ( spin[i] == VACANT ) zeros++;
		int ind_end = xyz_to_local(new_pos[0], new_pos[1], new_pos[2]);
		
//		if ( spin[ind_end] != FRAME ) {
//			printf("proc: %d ind_end: %d (global: %d) vac_new_pos: (%lf, %lf, %lf) grid: (%d, %d, %d) xyz[ind_end] (%lf, %lf, %lf) actual_state: %d\n", me, ind_end, id[ind_end], new_pos[0], new_pos[1], new_pos[2], gridx_end, gridy_end, gridz_end, xyz[ind_end][0], xyz[ind_end][1], xyz[ind_end][2], spin[ind_end]); 
//		}
		assert( spin[ind_end] == FRAME );	
	
		p[0] = (xyz[ind_end][0] - xyz[ind_vac][0]) / limit_t;
		p[1] = (xyz[ind_end][1] - xyz[ind_vac][1]) / limit_t;
		p[2] = (xyz[ind_end][2] - xyz[ind_vac][2]) / limit_t;
	
		int stepsx = abs(gridx_end - gridx_vac);
		int stepsy = abs(gridy_end - gridy_vac);
		int stepsz = abs(gridz_end - gridz_vac);
	
		int steps = MAX( stepsx, MAX(stepsy, stepsz) );
	
		if ( steps < 2 )  // Just one step means to exchange vacancy and frame -> NO ! it destroys frame
			return;
	
		double delta_tt = limit_t / steps;
	
		double x, y, z;
		int current_index=ind_vac, next_index;
		x = xyz[ind_vac][0];
		y = xyz[ind_vac][1];
		z = xyz[ind_vac][2];

		for ( int i = 0; i < steps; i++ ) {
			x += delta_tt * p[0];
			y += delta_tt * p[1];
			z += delta_tt * p[2];
			next_index = xyz_to_local(x, y, z);
			assert( next_index != current_index );
			if ( spin[next_index] > FRAME && next_index != ind_end ) {
				spin[current_index] = spin[next_index];
			}
			else {
				spin[current_index] = VACANT;
				break;
			}	
			current_index = next_index;	
		}
//		int after_zeros=0;
//		for ( int i = 0; i < nlocal; ++i )
//			if ( lattice[i] == VACANT ) after_zeros++;
//		assert ( zeros == after_zeros );
	}
	else {	
		// approximate xyz stored for final position
		double x_end = domain->boxxlo + gridx_end * Dx;
		double y_end = domain->boxylo + gridy_end * Dy;
		double z_end = domain->boxzlo + gridz_end * Dz;
	
		p[0] = (x_end - xyz[ind_vac][0]) / limit_t;
		p[1] = (y_end - xyz[ind_vac][1]) / limit_t;
		p[2] = (z_end - xyz[ind_vac][2]) / limit_t;
	
		int stepsx = abs(gridx_end - gridx_vac);
		int stepsy = abs(gridy_end - gridy_vac);
		int stepsz = abs(gridz_end - gridz_vac);
	
		int steps = MAX( stepsx, MAX(stepsy, stepsz) );
	
		if ( steps < 2 ) // Just one step means to exchange vacancy and frame -> NO ! it destroys frame
			return;
			if ( size_collapsinglist > (collapsing_directions.size() - 40 ) ) {
				try {
					collapsing_directions.resize( collapsing_directions.size() * 10 );
				}
				catch(std::bad_alloc& ba) {
					error->all(FLERR,"bad_alloc caught inside register_collapsing_event");
				}
			}
			collapsing_directions[size_collapsinglist++] = xyz[ind_vac][0]; // x
			collapsing_directions[size_collapsinglist++] = xyz[ind_vac][1]; // y
			collapsing_directions[size_collapsinglist++] = xyz[ind_vac][2]; // z
			collapsing_directions[size_collapsinglist++] = p[0]; // p.x		
			collapsing_directions[size_collapsinglist++] = p[1]; // p.y	
			collapsing_directions[size_collapsinglist++] = p[2]; // p.z	
			collapsing_directions[size_collapsinglist++] = limit_t; // Total t	
			collapsing_directions[size_collapsinglist++] = steps;
			collapsing_directions[size_collapsinglist++] = owner[ind_vac]; // Vacancy is owned by other processor
	}	
}

/* ----------------------------------------------------------------------
   In Parallel execution, consolidate the list of spins for the 
	annihilations to be made by each processor. Use the spin to calculate 
	the center of mass for each adjacent grain (grain with that spin).
	The average for the center of mass is not calculated, instead the
	sum of xyz positions and the number of sites per grain is consolidated
	and returned in vector cm_grains
------------------------------------------------------------------------- */
/*
void AppSinter::calculate_mass_center_distributed_grains( const int SIZE_LIST, const int SIZE_GROUP, vector<double> & cm_grains )
{ 
	const int SIZE_BUFFER ( SIZE_LIST*SIZE_GROUP );
	
	vector<int> annihilation_buffer( SIZE_LIST );
	vector<double> cm_buffer ( SIZE_BUFFER );
				
	MPI_Allgather(&size_annihilist, 1, MPI_INT, &annihilist_dist[0], 1, MPI_INT, world);
				
	int sum=0;
	for ( int i = 0; i < nprocs; i++ ) {
		displacement[i] = sum;
		sum += annihilist_dist[i];
	}
				
	MPI_Allgatherv(&annihilation_spin[0], size_annihilist, MPI_INT, &annihilation_buffer[0],
	&annihilist_dist[0], &displacement[0], MPI_INT, world);
			
	for ( int i = 0; i < SIZE_BUFFER; i++ ) {
		cm_buffer[i] = cm_grains[i] = 0.;
	}	

	for ( int i = 0; i < nlocal; i++ ) {
		for( int j = 0; j < SIZE_LIST; j++ ) {
			if ( spin[i] == annihilation_buffer[j] ) {
				for ( int k = 0; k < dimension; k++ )
					cm_buffer[j*SIZE_GROUP+k] += xyz[i][k];
				(cm_buffer[j*SIZE_GROUP+dimension])++;	
			}
		}
	}
	MPI_Allreduce(&cm_buffer[0], &cm_grains[0], SIZE_BUFFER, MPI_DOUBLE, MPI_SUM, world);
}
*/
/* ----------------------------------------------------------------------
   In Parallel execution, after the center of mass of the adjcent grains
	for all annihilation events is calculated, each processor compute the 
	vacancy new position and generate a list of collapsing paths. 
	The collapsing list is calculated by the processor that owns
	the vacancy to be annihilated, thus it is a local list.
	If collapsing some path moves one of the vacancies registered to be 
	annihilated, then this annihilation is removed from the list and no 
	collapsing is computed.
	The list of local pending collapsing paths is stored in vector:
	collapsing_directions
------------------------------------------------------------------------- */

/*
void AppSinter::generate_local_vacancy_collapsing_list( const int SIZE_LIST, const int SIZE_GROUP, vector<double> & cm_grains )		
{
	const int SIZE_BUFFER ( SIZE_LIST*SIZE_GROUP );
	
	vector<double> new_pos ( dimension );
	vector<double> p ( dimension );
	double mint;
	
	int local_index = displacement[me] * SIZE_GROUP;
	vector<double> coord_cm ( dimension );
	
	for ( int i = 0; i < size_annihilist; i++ ) {
		int j = annihilation_list[i];
		if ( spin[j] != VACANT ) { // Other vacancy update modified this site -> no longer a vacancy
			naccept--;
			vac_made--;
			local_index += SIZE_GROUP;
			continue;
		}	
		
		assert ( local_index >= 0 && local_index < SIZE_BUFFER );
		for ( int k = 0; k < dimension; k++ )
			coord_cm[k] = cm_grains[local_index+k] / cm_grains[local_index+dimension];
			
			calculate_vacancy_new_position( j, coord_cm, p, mint, new_pos );
			register_collapsing_event(j, p, mint, new_pos);
			local_index += SIZE_GROUP;
			}
}
*/

void AppSinter::generate_local_vacancy_collapsing_list( const int SIZE_LIST )		
{
	vector<double> coord_cm ( dimension );
	vector<double> new_pos ( dimension );
	vector<double> p ( dimension );
	double mint;
	
	for ( int i = 0; i < size_annihilist; i++ ) {
		int j = annihilation_list[i];
		if ( spin[j] != VACANT ) { // Other vacancy update modified this site -> no longer a vacancy
			naccept--;
			vac_made--;
			continue;
		}	
		
		// Find the spin of the adjacent grain in the global list of grains
		// and retrieve the center of mass
		int grain_ind = -1;
		for ( int ind = 0; ind < numgrains; ind++ ) {
			if ( annihilation_spin[i] == grain_spins[ind] ) {
				grain_ind = ind;
				break;
			}
		}
		/*		if ( grain_ind == -1 ) {
		 printf("proc: %d vacancy to be annihilated: %d (spin: %d) spin of adjacent grain: %d\n", me, j, lattice[j], annihilation_spin[i]);
		 printf("List of neighboring spins:\n");
		 for ( int count = 0; count < numneigh[j]; count++ ) {
		 printf("%d ", lattice[neighbor[j][count]] );
		 }
		 printf("\nnumber of grains: %d. List of spins:\n", numgrains );  
		 //			for ( int count = 0; count < numgrains; count++ ) {
		 //				printf("%d ", grain_spins[count]);
		 //			}
		 //			printf("\n");
		 }
		 */		
		assert ( grain_ind > -1 ); // The adjacent grain was found !
		
		int offset = grain_ind * dimension;		
		for ( int k = 0; k < dimension; k++ )
			coord_cm[k] = grain_mass_center[offset+k];
		
		calculate_vacancy_new_position( j, coord_cm, p, mint, new_pos );
		register_collapsing_event(j, p, mint, new_pos);
	}
}			

/* ----------------------------------------------------------------------
   In Parallel execution, after some collapsing paths are calculated
	locally it could happen that some vacancies are shifted. Thus, the
	local collapsed list is reviewed and vacancies to be annihilated that
	are not longer there are removed from the list.
	The list of local pending collapsing paths is stored in vector:
	collapsing_directions
------------------------------------------------------------------------- */

void AppSinter::filter_local_vacancy_collapsing_list()
{
	vector<double> copy_collapsing_directions ( size_collapsinglist );
		
	for ( int i = 0; i < size_collapsinglist; i++ )
		copy_collapsing_directions[i] = collapsing_directions[i];

	int m = 0, j = 0, ilocalt, stepst, ownert;	
	double xt, yt, zt, pxt, pyt, pzt, limittt;
	bool copyflag;
	while ( m < size_collapsinglist ) {
		xt = copy_collapsing_directions[m++];
		yt = copy_collapsing_directions[m++];
		zt = copy_collapsing_directions[m++];
		pxt = copy_collapsing_directions[m++];
		pyt = copy_collapsing_directions[m++];
		pzt = copy_collapsing_directions[m++];
		limittt = copy_collapsing_directions[m++];
		stepst = (int)(copy_collapsing_directions[m++]);
		ownert = (int)(copy_collapsing_directions[m++]);
		copyflag = false;
		if ( ownert != me ) {
			copyflag = true;
		}	
		else {
			ilocalt = xyz_to_local(xt, yt, zt);
			if ( spin[ilocalt] == VACANT ) 
				copyflag = true;
		}
		if ( copyflag ) {		
			collapsing_directions[j++] = xt;
			collapsing_directions[j++] = yt;
			collapsing_directions[j++] = zt;
			collapsing_directions[j++] = pxt;
			collapsing_directions[j++] = pyt;
			collapsing_directions[j++] = pzt;
			collapsing_directions[j++] = limittt;
			collapsing_directions[j++] = stepst;
		}
	}		
	size_collapsinglist = j;
}

/* ----------------------------------------------------------------------
   In Parallel execution, consolidate the list of pending annihilations,
	follow the different paths updating the sites and exchange boundary 
	information.
	All the processors folllow the paths but just the one who ownes the 
	sites makes the corresponding updates. When a border between processors 
	is found, the previous processor keeps waiting for the new state comming
	from the next processor. The information exchanged is global index and
	spin update.
	If the boder (frame) of the simulation space is found before completing
	the supposed number of steps and the expected final site is in a 
	different processor, then the processor waits for some information (to
	avoid blocking) but does not use it, as it is just frame sites.
	The border states are collected in a stack that is updated at the
	end of each path.
	If a global index corresponding to a vacancy to be annihilated has been 
	updated (as indicated by a list of iglobal indices updated every time
	a new site is sweeped), then the corresponding collapsing route is not
	comlpeted.
	To avoid destroying the border, before updating any site it is veryfied
	that it does not correspond to the border. 
------------------------------------------------------------------------- */

void AppSinter::collapse_pending_list ( const int SIZE_LIST ) 
{
	vector<double> collapsing_directions_buffer ( SIZE_LIST );
			
	vector <int> collapsing_dist ( nprocs );
	MPI_Allgather(&size_collapsinglist,1,MPI_INT,&collapsing_dist[0],1,MPI_INT,world);
				
	int sum=0;
	for ( int i = 0; i < nprocs; i++ ) {
		displacement[i] = sum;
		sum += collapsing_dist[i];
	}
	MPI_Allgatherv(&collapsing_directions[0], size_collapsinglist, MPI_DOUBLE,
	&collapsing_directions_buffer[0], &collapsing_dist[0], &displacement[0], MPI_DOUBLE, world);

	int steps, proc, next_proc, new_state, iglobal, next_iglobal, next_index, current_index;
	double x0, y0, z0, px, py, pz, limit_t, x, y, z;
	std::map<int,int>::iterator loc;
				
	std::stack<int> pending_updates;	
				
	std::vector<int> included;
								
	int border_send[2], border_recv[2];
	MPI_Status status;
				
	// Buffer chunks: from (x,y,z), direction (px, py, pz), limit_t, steps
	int m = 0;
	while ( m < SIZE_LIST ) {
		x0 = static_cast<double> (collapsing_directions_buffer[m++]);
		y0 = static_cast<double> (collapsing_directions_buffer[m++]);
		z0 = static_cast<double> (collapsing_directions_buffer[m++]);
		px = static_cast<double> (collapsing_directions_buffer[m++]);
		py = static_cast<double> (collapsing_directions_buffer[m++]);
		pz = static_cast<double> (collapsing_directions_buffer[m++]);
		limit_t = static_cast<double> (collapsing_directions_buffer[m++]);
		steps = (int)(static_cast<double> (collapsing_directions_buffer[m++]));
					
		x = x0;
		y = y0;
		z = z0;

		bool duplicate = false;
		proc = xyz_to_processor_and_global( x, y, z, iglobal );
		for ( int i = 0; i < included.size(); i++ ) {
			if ( included[i] == iglobal ) {
				duplicate = true;
			}	
		}
		if ( !duplicate )		
			included.push_back( iglobal );
		else {
			if ( proc == me ) {
				naccept--;
				vac_made--;
			}
			continue;	
		}
		
		int dum_ig, dum_local_i=-1, dum_state;
		if ( proc == me ) {
			loc = hash.find( iglobal );
			current_index = loc->second;
			assert( spin[current_index] == VACANT );
//			if ( spin[current_index] > VACANT ) {
//				printf("proc: %d collapse start in iglobal: %d ilocal: %d state: %d\n", me, iglobal, current_index, spin[current_index]); 
//			}							
		}
		else
			current_index = -1;
					
		double delta_tt = limit_t / steps;
		bool breakflag = 0;
		int last_proc_breakflag, next_proc_breakflag;
				
		for ( int i = 0; i < steps; i++ ) {
			x += delta_tt * px;
			y += delta_tt * py;
			z += delta_tt * pz;
			next_proc = xyz_to_processor_and_global( x, y, z, next_iglobal );
			included.push_back( next_iglobal );
			if ( next_proc == me ) {
				loc = hash.find(next_iglobal);
				next_index = loc->second;
				if ( current_index != -1 ) { // Update inside same processor
					if ( spin[next_index] > FRAME ) {
						spin[current_index] = spin[next_index];
					}
					else {
						spin[current_index] = VACANT;
						int dum;
						last_proc_breakflag = xyz_to_processor_and_global(x0+limit_t*px, y0+limit_t*py, z0+limit_t*pz, dum);
						if ( last_proc_breakflag != me ) {
							do {
								x += delta_tt * px;
								y += delta_tt * py;
								z += delta_tt * pz;
								next_proc = xyz_to_processor_and_global( x, y, z, dum );
							} while ( next_proc == me );
							next_proc_breakflag = next_proc;
							breakflag = 1;
						}		
						break;
					}
				}
				else { // Entering a new processor
					border_send[0] = iglobal;
					if ( spin[next_index] > FRAME ) {
						border_send[1] = spin[next_index];
					}
					else {
						border_send[1] = VACANT;
						int dum;
						last_proc_breakflag = xyz_to_processor_and_global(x0+limit_t*px, y0+limit_t*py, z0+limit_t*pz, dum);
						if ( last_proc_breakflag != me ) {
							do {
								x += delta_tt * px;
								y += delta_tt * py;
								z += delta_tt * pz;
								next_proc = xyz_to_processor_and_global( x, y, z, dum );
							} while ( next_proc == me );
							next_proc_breakflag = next_proc;
							breakflag = 1;
						}		
						i = steps; // force end of for after sending update
					}	
					MPI_Send(&border_send[0], 2, MPI_INT, proc, 0, world);
				}
				current_index = next_index;	
			}
			else if ( proc == me ) { // Leaving this processor
				MPI_Recv(&border_recv[0], 2, MPI_INT, next_proc, 0, world, &status);
				pending_updates.push(border_recv[1]); // new_state
				pending_updates.push(border_recv[0]); // iglobal
			}
			proc = next_proc;
			iglobal = next_iglobal;
		}
		if ( breakflag ) {			
			MPI_Recv(&border_recv[0], 2, MPI_INT, next_proc_breakflag, 0, world, &status);						
		}
		int ilocal, new_state;
					
		while ( pending_updates.size() ) {
			iglobal = pending_updates.top();
			pending_updates.pop();
			new_state = pending_updates.top();
			pending_updates.pop();											
			loc = hash.find(iglobal);
			ilocal = loc->second;					
			if ( spin[ilocal] != FRAME )
				spin[ilocal] = new_state;
		}				
	}
	
	size_collapsinglist = 0;															
}

/////////////////////////////////////////// GRAIN BOUNDARY CALCULATIONS /////////////////////////////////////////////////	

double AppSinter::calculate_gb_update(double current_time)
{
	double grain_size_avg;
	double prev_density = -1;
	
	if ( gb_ini < 0.0 ) {
		gb_ini = calculate_gb_average(grain_size_avg);
		gsize_ini = grain_size_avg;
		delta_t = 20;
//		if( me == 0 ) printf("grain size average at beginning of sintering: %lf\n",grain_size_avg);
	}
	else {
		double gb_update = calculate_gb_average(grain_size_avg);
		gb_factor = (gb_ini * gb_ini) / (gb_update * gb_update);
		prev_density = density; 
//		if ( me == 0 ) {
//			printf("gb_factor: %lf\n", gb_factor);
//			printf("grain size average: %lf\n",grain_size_avg);
//		}
//		if ( me == 0 )
//			printf("density: %lf\n", density );
	}
	
  density = calculate_density();
	
	
//	double sintertime =	stoptime - time_sinter_start;
//	double currentoffset = current_time - time_sinter_start;
	double delta_density = fabs(density - prev_density);
	
	if ( delta_density > 0.05 )
	  delta_t /= 2;
	else if ( delta_density < 0.01 )
	  delta_t *= 2;
/*	
	if ( currentoffset < 0.05 * sintertime )
		delta_t = 0.01 * sintertime;
	else if ( currentoffset < 0.1 * sintertime )
		delta_t = 0.025 * sintertime;
	else if ( currentoffset < 0.35 * sintertime )
		delta_t = 0.05 * sintertime;
	else
		delta_t = 0.15 * sintertime;
*/
/*	
	if ( currentoffset < 0.05 * sintertime )
		delta_t = 0.01 * sintertime;
	else if ( currentoffset < 0.1 * sintertime )
		delta_t = 0.025 * sintertime;
	else if ( currentoffset < 0.35 * sintertime )
		delta_t = 0.05 * sintertime;
	else
		delta_t = 0.10 * sintertime;
*/	
	return ( current_time + delta_t );		
}

/* ----------------------------------------------------------------------
   Calculate grain boundary average to compute frequency of annihilations
------------------------------------------------------------------------- */
/*
double AppSinter::calculate_gb_average(double & grain_size_average)
{  // collects all in a central processor, per grain, for final calc.
	int localsum_grain_sites = 0;
	int localsum_face_sites = 0;
	int localsum_faces = 0;
	int localnum_grains = 0;
	
	vector<bool> site_included( nlocal+nghost, false );
	
	vector<int> grains_spins;
	int num_faces=0;
	
	vector<int> local_buffer( nlocal * 10 );
	int size_local_buffer = 0;
	
	int num_grains = 0;
	for ( int i = 0; i < nlocal; i++ ) {
		if ( spin[i] == VACANT || spin[i] == FRAME )
			continue;
		if ( site_included[i] )
			continue;
		grains_spins.push_back( spin[i] );	
		vector<int> faces;
		vector<int> neigh_procs;
		int grain_vol, face_sites;
		bool multiproc;
		
		cluster_faces( i, site_included, grain_vol, face_sites, faces, multiproc, neigh_procs );
		if ( multiproc ) {
			local_buffer[size_local_buffer++] = spin[i];
			local_buffer[size_local_buffer++] = faces.size();
			for ( int j = 0; j < faces.size(); j++ ) {
				local_buffer[size_local_buffer++] = faces[j];
			}
			local_buffer[size_local_buffer++] = neigh_procs.size();
			for ( int j = 0; j < neigh_procs.size(); j++ ) {
				local_buffer[size_local_buffer++] = neigh_procs[j];
			}
		}	
		else {
			localsum_faces += faces.size();
			localnum_grains++;
		}
		localsum_face_sites += face_sites;
		localsum_grain_sites += grain_vol;	
	}
	local_buffer[size_local_buffer++] = localsum_faces;
	local_buffer[size_local_buffer++] = localsum_face_sites;
	local_buffer[size_local_buffer++] = localsum_grain_sites;
	local_buffer[size_local_buffer++] = localnum_grains;
	
	double gb_avg = 0;
	grain_size_average = 0;

	// Consolidate buffer in proc 0
	int size_buffer_all = 0, tmp=0, nrecv, m, tmp_all;
	MPI_Reduce(&size_local_buffer,&size_buffer_all,1,MPI_INT,MPI_SUM,0,world);
	MPI_Status status;
	MPI_Request request;
	
	int total_faces = 0, total_face_sites = 0, total_grain_sites = 0, total_grains = 0;
	
	if ( me == 0 ) {
		vector<int> copy_buffer( size_buffer_all );
		size_local_buffer -= 4;
		int offset = size_local_buffer;	
		total_faces = localsum_faces;
		total_face_sites = localsum_face_sites;
		total_grain_sites = localsum_grain_sites;
		total_grains = localnum_grains;
		
		vector<int> rec_data( nprocs, 0 );
		for (int iproc = 1; iproc < nprocs; iproc++) {
			MPI_Irecv(&copy_buffer[offset], size_buffer_all-offset, MPI_INT, iproc, 0, world, &request);
			MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
			MPI_Wait(&request,&status);
			MPI_Get_count(&status, MPI_INT, &nrecv);
			total_faces += copy_buffer[offset+nrecv-4];
			total_face_sites += copy_buffer[offset+nrecv-3];
			total_grain_sites += copy_buffer[offset+nrecv-2];
			total_grains += copy_buffer[offset+nrecv-1];
			rec_data[iproc] = nrecv - 4 + rec_data[iproc-1];
			offset += nrecv - 4;
		}	
		// copy local buffer
		for ( int i = 0; i < size_local_buffer; i++ ) {
			copy_buffer[i] = local_buffer[i];
		}	
		
		int spin, num_faces, neigh_procs, next_proc;
		m = 0;
		while (m < offset) {
			int ind_current_spin = m;
			spin = copy_buffer[m++];
			num_faces = copy_buffer[m++];
			vector<int> faces( num_faces );
			for ( int i = 0; i < num_faces; i++ ) {
				faces[i] = copy_buffer[m++];
			}
			neigh_procs = copy_buffer[m++];
			for ( int i = 0; i < neigh_procs; i++ ) {
				next_proc = copy_buffer[m++];
			}
			if ( spin == -1 ) // grain burned
				continue;
			int k = m, aux_faces, aux_neighproc;
			while ( k < offset ) {
				int ospin = copy_buffer[k++];
				if ( ospin != spin ) {
					aux_faces = copy_buffer[k];
					k += (aux_faces + 1);
					aux_neighproc = copy_buffer[k];
					k += (aux_neighproc + 1);
					continue;
				}							
				copy_buffer[k-1] = -1; // burn grain in buffer
				aux_faces = copy_buffer[k++];
				for ( int i = 0; i < aux_faces; i++ ) {
					faces.push_back( copy_buffer[k++] );
				}	
				aux_neighproc = copy_buffer[k];
				k += (aux_neighproc + 1);
			}
			// Count number of faces by removing duplicates from list of faces and burning face already included 
			int count_faces = 0;
			for ( int i = 0; i < faces.size(); i++ ) {
				if ( faces[i] != -1 ) {
					int neigh_face = faces[i];
					count_faces++;
					faces[i] = -1;
					for ( int j = i+1; j < faces.size(); j++ ) {
						if ( faces[j] == neigh_face ) {
							faces[j] = -1;
						}	
					}	
				}	
			}
			total_faces += count_faces;
			total_grains++;
			copy_buffer[ind_current_spin] = -1; // burning current spin 
		}	
		gb_avg = (double) total_face_sites / (double) total_faces;
		grain_size_average = (double) total_grain_sites / (double) total_grains;
	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
		MPI_Rsend(&local_buffer[0],size_local_buffer,MPI_INT,0,0,world);
	}
	
	vector<double> avg(2);
	if ( me == 0 ) {
		avg[0] = gb_avg;
		avg[1] = grain_size_average;
	}
	 
	MPI_Bcast(&avg[0], 2, MPI_DOUBLE, 0, world);

	if ( me > 0 ) {
		gb_avg = avg[0];
		grain_size_average = avg[1];
//		printf("proc: %d gb_avg: %lf grain_size_average: %lf\n", me, gb_avg, grain_size_average ); 		
	}
	return gb_avg;
}
*/

double AppSinter::calculate_gb_average(double & grain_size_average)
{
	int localsum_grain_sites = 0;
	int localsum_face_sites = 0;
	int localsum_faces = 0;
	int localnum_grains = 0;
	
	vector<bool> site_included( nlocal+nghost, false );
	
	vector<int> grains_spins;
	int num_faces=0;
	
	vector<int> local_buffer( nlocal * 10 );
	int size_local_buffer = 0;
	
	const int GROUPSIZE ( dimension + 1 );
	vector<double> center_mass( dimension );
	vector<int> local_grain_spin;
	vector<double> local_grain_mass_center;
	int size_local_numgrain = 0;
	if ( nprocs > 1 ) {
		local_grain_spin.resize( nlocal );
		local_grain_mass_center.resize( nlocal * GROUPSIZE );
	}
	
	int num_grains = 0;
	for ( int i = 0; i < nlocal; i++ ) {
		if ( spin[i] == VACANT || spin[i] == FRAME )
			continue;
		if ( site_included[i] )
			continue;
		grains_spins.push_back( spin[i] );	
		vector<int> faces;
		vector<int> neigh_procs;
		int grain_vol, face_sites;
		bool multiproc;
		
		cluster_faces( i, site_included, grain_vol, face_sites, faces, multiproc, neigh_procs, center_mass );
		if ( multiproc ) {
			local_buffer[size_local_buffer++] = spin[i];
			local_buffer[size_local_buffer++] = faces.size();
			for ( int j = 0; j < faces.size(); j++ ) {
				local_buffer[size_local_buffer++] = faces[j];
			}
			local_buffer[size_local_buffer++] = neigh_procs.size();
			for ( int j = 0; j < neigh_procs.size(); j++ ) {
				local_buffer[size_local_buffer++] = neigh_procs[j];
			}
		}	
		else {
			localsum_faces += faces.size();
			localnum_grains++;
		}
		localsum_face_sites += face_sites;
		localsum_grain_sites += grain_vol;
		
		// Store spin and center of mass of each grain
		if ( nprocs > 1 ) {
			local_grain_spin[size_local_numgrain] = spin[i];
			for ( int k = 0; k < dimension; k++ )
				local_grain_mass_center[size_local_numgrain*GROUPSIZE+k] = center_mass[k];
			local_grain_mass_center[size_local_numgrain*GROUPSIZE+dimension] = grain_vol;
			size_local_numgrain++;
		}	
	}
	local_buffer[size_local_buffer++] = localsum_faces;
	local_buffer[size_local_buffer++] = localsum_face_sites;
	local_buffer[size_local_buffer++] = localsum_grain_sites;
	local_buffer[size_local_buffer++] = localnum_grains;
	
	double gb_avg = 0;
	grain_size_average = 0;
	
	// Consolidate buffer in proc 0
	int size_buffer_all = 0, tmp=0, nrecv, m, tmp_all;
	MPI_Reduce(&size_local_buffer,&size_buffer_all,1,MPI_INT,MPI_SUM,0,world);
	MPI_Status status;
	MPI_Request request;
	
	int total_faces = 0, total_face_sites = 0, total_grain_sites = 0, total_grains = 0;
	
	// Consolidate list of grains mass center
	if ( nprocs > 1 ) {
		//		if ( me == 0 ) printf("going to calculate center of mass at time: %lf\n", time );
		consolidate_mass_center_distributed_grains( local_grain_spin, local_grain_mass_center, size_local_numgrain );	
		//		printf("after calculating center of mass. proc: %d local number of grains: %d total: %d\n", me, size_local_numgrain, numgrains);
	}
	
	// Calculations regarding grain boundary
	if ( me == 0 ) {
		vector<int> copy_buffer( size_buffer_all );
		size_local_buffer -= 4;
		int offset = size_local_buffer;	
		total_faces = localsum_faces;
		total_face_sites = localsum_face_sites;
		total_grain_sites = localsum_grain_sites;
		total_grains = localnum_grains;
		
		vector<int> rec_data( nprocs, 0 );
		for (int iproc = 1; iproc < nprocs; iproc++) {
			MPI_Irecv(&copy_buffer[offset], size_buffer_all-offset, MPI_INT, iproc, 0, world, &request);
			MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
			MPI_Wait(&request,&status);
			MPI_Get_count(&status, MPI_INT, &nrecv);
			total_faces += copy_buffer[offset+nrecv-4];
			total_face_sites += copy_buffer[offset+nrecv-3];
			total_grain_sites += copy_buffer[offset+nrecv-2];
			total_grains += copy_buffer[offset+nrecv-1];
			rec_data[iproc] = nrecv - 4 + rec_data[iproc-1];
			offset += nrecv - 4;
		}	
		// copy local buffer
		for ( int i = 0; i < size_local_buffer; i++ ) {
			copy_buffer[i] = local_buffer[i];
		}	
		
		int spin, num_faces, neigh_procs, next_proc;
		m = 0;
		while (m < offset) {
			int ind_current_spin = m;
			spin = copy_buffer[m++];
			num_faces = copy_buffer[m++];
			vector<int> faces( num_faces );
			for ( int i = 0; i < num_faces; i++ ) {
				faces[i] = copy_buffer[m++];
			}
			neigh_procs = copy_buffer[m++];
			for ( int i = 0; i < neigh_procs; i++ ) {
				next_proc = copy_buffer[m++];
			}
			if ( spin == -1 ) // grain burned
				continue;
			int k = m, aux_faces, aux_neighproc;
			while ( k < offset ) {
				int ospin = copy_buffer[k++];
				if ( ospin != spin ) {
					aux_faces = copy_buffer[k];
					k += (aux_faces + 1);
					aux_neighproc = copy_buffer[k];
					k += (aux_neighproc + 1);
					continue;
				}							
				copy_buffer[k-1] = -1; // burn grain in buffer
				aux_faces = copy_buffer[k++];
				for ( int i = 0; i < aux_faces; i++ ) {
					faces.push_back( copy_buffer[k++] );
				}	
				aux_neighproc = copy_buffer[k];
				k += (aux_neighproc + 1);
			}
			// Count number of faces by removing duplicates from list of faces and burning face already included 
			int count_faces = 0;
			for ( int i = 0; i < faces.size(); i++ ) {
				if ( faces[i] != -1 ) {
					int neigh_face = faces[i];
					count_faces++;
					faces[i] = -1;
					for ( int j = i+1; j < faces.size(); j++ ) {
						if ( faces[j] == neigh_face ) {
							faces[j] = -1;
						}	
					}	
				}	
			}
			total_faces += count_faces;
			total_grains++;
			copy_buffer[ind_current_spin] = -1; // burning current spin 
		}	
		gb_avg = (double) total_face_sites / (double) total_faces;
		grain_size_average = (double) total_grain_sites / (double) total_grains;
	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
		MPI_Rsend(&local_buffer[0],size_local_buffer,MPI_INT,0,0,world);
	}
	
	vector<double> avg(2);
	if ( me == 0 ) {
		avg[0] = gb_avg;
		avg[1] = grain_size_average;
	}
	
	MPI_Bcast(&avg[0], 2, MPI_DOUBLE, 0, world);
	
	if ( me > 0 ) {
		gb_avg = avg[0];
		grain_size_average = avg[1];
		//		printf("proc: %d gb_avg: %lf grain_size_average: %lf\n", me, gb_avg, grain_size_average ); 		
	}
	return gb_avg;
}


/* ---------------------------------------------------------------------- */
/*
void AppSinter::cluster_faces( int start_ilocal, vector<bool> & site_included, int & grain_vol, int & face_sites, vector<int> & faces,
										 bool & multiproc, vector<int> & neigh_procs )
{  // per grain per processor
	int lspin = spin[start_ilocal];

	face_sites = 0;
	multiproc = false;
	
	std::stack<int> exploring;
	
	exploring.push( start_ilocal );
	site_included[start_ilocal] = true;
	grain_vol = 1;
	
	vector<int> aux_faces;
	vector<int> aux_neigh_procs; // To store neighboring processors that store part of the grain
	
	while ( exploring.size() ) {
		int ilocal = exploring.top();
		exploring.pop();
		int nevent = 0, m;
		for ( int nbor = 0; nbor < numneigh[ilocal]; nbor++ ) {
			int neigh = neighbor[ilocal][nbor];
			int neigh_spin = spin[neigh];
			if ( neigh_spin == FRAME || neigh_spin == VACANT ) continue; // Outside simulation space or vacancy
			if ( neigh_spin == lspin ) { // Same spin, part of the cluster 
				if ( !site_included[neigh] ) { // if not included --> include in exploring stack
					exploring.push( neigh );
					site_included[neigh] = true;
					if ( neigh < nlocal )	grain_vol++;
					else {
						aux_neigh_procs.push_back ( owner[neigh] );
						multiproc = true; // grain is contained in more than one processor
					}	
				}	
				continue;
			}	
			if ( ilocal >= nlocal ) continue; // It is going to be counted in the processor that owns it
			for ( m = 0; m < nevent; m++ )
				if ( neigh_spin == unique[m] ) break; // Registered Event
			if ( m < nevent ) continue; // Previous cycle was interrupted because event was found
			unique[nevent++] = neigh_spin;
		}
		face_sites += nevent;
		for ( int i = 0; i < nevent; i++ ) {
				aux_faces.push_back( unique[i] );
		}	
	}
	// Removing duplicates from list of faces by burning face already included 
	for ( int i = 0; i < aux_faces.size(); i++ ) {
		if ( aux_faces[i] != -1 ) {
			int neigh_spin = aux_faces[i];
			faces.push_back( neigh_spin );
			aux_faces[i] = -1;
			for ( int j = i+1; j < aux_faces.size(); j++ ) {
				if ( aux_faces[j] == neigh_spin ) {
					aux_faces[j] = -1;
				}	
			}	
		}	
	}

	if ( multiproc ) { // Grain is contained in more than one processor
		// Removing duplicates from list of neighboring procs by burning proc already included 
		for ( int i = 0; i < aux_neigh_procs.size(); i++ ) {
			if ( aux_neigh_procs[i] != -1 ) {
				int neighproc = aux_neigh_procs[i];
				neigh_procs.push_back( neighproc );
				aux_neigh_procs[i] = -1;
				for ( int j = i+1; j < aux_neigh_procs.size(); j++ ) {
					if ( aux_neigh_procs[j] == neighproc ) {
						aux_neigh_procs[j] = -1;
					}	
				}	
			}	
		}
	}
} 
*/

/* ----------------------------------------------------------------------
 register per grain: sites in faces, in edges and in corners.
 For sites in faces register spin of nieghboring grain, thus if
 grain is distributed in several processors, just different faces are
 counted.
 ------------------------------------------------------------------------- */

void AppSinter::cluster_faces( int start_ilocal, vector<bool> & site_included, int & grain_vol, int & face_sites, vector<int> & faces,
							  bool & multiproc, vector<int> & neigh_procs, vector<double> & cm )
{
	int ispin = spin[start_ilocal];
	
	face_sites = 0;
	multiproc = false;
	
	std::stack<int> exploring;
	
	for ( int k = 0; k < dimension; k++ )
		cm[k] = 0.0;
	
	exploring.push( start_ilocal );
	site_included[start_ilocal] = true;
	grain_vol = 1;
	
	vector<int> aux_faces;
	vector<int> aux_neigh_procs; // To store neighboring processors that store part of the grain
	
	while ( exploring.size() ) {
		int ilocal = exploring.top();
		exploring.pop();
		int nevent = 0, m;
		for ( int nbor = 0; nbor < numneigh[ilocal]; nbor++ ) {
			int neigh = neighbor[ilocal][nbor];
			int neigh_spin = spin[neigh];
			if ( neigh_spin == FRAME || neigh_spin == VACANT ) continue; // Outside simulation space or vacancy
			if ( neigh_spin == ispin ) { // Same spin, part of the cluster 
				if ( !site_included[neigh] ) { // if not included --> include in exploring stack
					exploring.push( neigh );
					site_included[neigh] = true;
					if ( neigh < nlocal )	{
						grain_vol++;
						for ( int k = 0; k < dimension; k++ )
							cm[k] += xyz[neigh][k];
					}
					else {
						aux_neigh_procs.push_back ( owner[neigh] );
						multiproc = true; // grain is contained in more than one processor
					}	
				}	
				continue;
			}	
			if ( ilocal >= nlocal ) continue; // It is going to be counted in the processor that owns it
			for ( m = 0; m < nevent; m++ )
				if ( neigh_spin == unique[m] ) break; // Registered Event
			if ( m < nevent ) continue; // Previous cycle was interrupted because event was found
			unique[nevent++] = neigh_spin;
		}
		face_sites += nevent;
		for ( int i = 0; i < nevent; i++ ) {
			aux_faces.push_back( unique[i] );
		}	
	}
	// Removing duplicates from list of faces by burning face already included 
	for ( int i = 0; i < aux_faces.size(); i++ ) {
		if ( aux_faces[i] != -1 ) {
			int neigh_spin = aux_faces[i];
			faces.push_back( neigh_spin );
			aux_faces[i] = -1;
			for ( int j = i+1; j < aux_faces.size(); j++ ) {
				if ( aux_faces[j] == neigh_spin ) {
					aux_faces[j] = -1;
				}	
			}	
		}	
	}
	
	if ( multiproc ) { // Grain is contained in more than one processor
		// Removing duplicates from list of neighboring procs by burning proc already included 
		for ( int i = 0; i < aux_neigh_procs.size(); i++ ) {
			if ( aux_neigh_procs[i] != -1 ) {
				int neighproc = aux_neigh_procs[i];
				neigh_procs.push_back( neighproc );
				aux_neigh_procs[i] = -1;
				for ( int j = i+1; j < aux_neigh_procs.size(); j++ ) {
					if ( aux_neigh_procs[j] == neighproc ) {
						aux_neigh_procs[j] = -1;
					}	
				}	
			}	
		}
	}
} 

/* ----------------------------------------------------------------------
 In Parallel execution, consolidate the list of centers of mass for
 each grain in the domain. The same spin comming from different
 processors is collected as one grain. 
 At the end every processor knows the center of mass of every grain.
 Spins are stored in vector 'grain_spins', center of mas in vector
 'grain_mass_center' and number of grain in 'numgrains'
 ------------------------------------------------------------------------- */

void AppSinter::consolidate_mass_center_distributed_grains( vector<int> & spins, vector<double> & cm_list, const int LOCALNUMGRAINS ) 
{
	const int GROUPSIZE ( dimension + 1 ); // coordinates plus one element for volume
	int local_size = LOCALNUMGRAINS; 
	int fat_numgrains=0;
	
	vector<int> cm_distribution( nprocs );
	
	MPI_Gather(&local_size, 1, MPI_INT, &cm_distribution[0], 1, MPI_INT, 0, world);
	
	vector<int> cm_displacement;
	vector<int> gspins;
	vector<double> gcm; 
	
	if ( me == 0 ) {			
		cm_displacement.resize( nprocs );
		int sum=0;
		for ( int i = 0; i < nprocs; i++ ) {
			cm_displacement[i] = sum;
			sum += cm_distribution[i];
		}
		
		gspins.resize( sum );
		gcm.resize( sum * GROUPSIZE ); 
		fat_numgrains = sum;
	}
	
	MPI_Gatherv(&spins[0], local_size, MPI_INT, &gspins[0], &cm_distribution[0], &cm_displacement[0], MPI_INT, 0, world);
	if ( me == 0 ) {
		for ( int i = 0; i < nprocs; i++ ) {
			cm_distribution[i] *= GROUPSIZE;
			cm_displacement[i] *= GROUPSIZE;
		}
	}
	MPI_Gatherv(&cm_list[0], GROUPSIZE*local_size, MPI_DOUBLE, &gcm[0], &cm_distribution[0], &cm_displacement[0], MPI_DOUBLE, 0, world);
	
	numgrains = 0;
	if ( me == 0 ) {
		grain_spins.resize( 0 );
		grain_mass_center.resize( 0 );
		for ( int i = 0; i < fat_numgrains; i++ ) {
			if ( gspins[i] != -1 ) { // Grain has not been counted
				grain_spins.push_back( gspins[i] );
				for ( int k = 0; k < dimension; k++ ) {
					grain_mass_center.push_back( gcm[i*GROUPSIZE+k] );
				}
				int gvol = gcm[i*GROUPSIZE+3];
				for ( int j = i+1; j < fat_numgrains; j++ ) {
					if ( gspins[j] == gspins[i] ) { // Another part of grain found
						for ( int k = 0; k < dimension; k++ ) {
							grain_mass_center[numgrains*dimension+k] += gcm[j*GROUPSIZE+k];
						}
						gvol += gcm[j*GROUPSIZE+3];
						gspins[j] = -1; // Burn grain because it has been counted
					}
				}
				// Calculate average
				for ( int k = 0; k < dimension; k++ ) {
					grain_mass_center[numgrains*dimension+k] /= gvol;
				}
				numgrains++;
				gspins[i] = -1; // Burn grain because it has been counted
			}
		}
	}
	
	MPI_Bcast( &numgrains, 1, MPI_INT, 0, world );
	if ( me > 0 ) {
		grain_spins.resize( numgrains );
		grain_mass_center.resize( numgrains * dimension );
	}
	MPI_Bcast( &grain_spins[0], numgrains, MPI_INT, 0, world );
	MPI_Bcast( &grain_mass_center[0], numgrains*dimension, MPI_DOUBLE, 0, world );	
	/*	
	 int i, j,count=0;
	 for ( i = 0; i < nlocal; i++ ) {
	 if ( lattice[i] > VACANT ) { 
	 for ( j = 0; j < numgrains; j++ ) {
	 if ( lattice[i] == grain_spins[j] )
	 break;
	 }
	 if ( lattice[i] != grain_spins[j] ) {
	 //	printf("proc: %d spin: %d not found in list of grains (j: %d)\n", me, lattice[i], j );
	 count++;
	 }
	 }	
	 }
	 printf("proc: %d number of grains not found in the list: %d\n", me, count );
	 printf("proc: %d size of vector grain_spins: %d\n", me, grain_spins.size() ); 
	 */	
}

/////////////////////////////////////////// STATISTICS /////////////////////////////////////////////////


/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppSinter::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %12s %14s %10s %12s","Time","Naccept","Nreject","Nsweeps", "Vmade");
}


/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppSinter::stats(char *strtmp)
{
  double naccept_double = (double)naccept;
  double naccept_double_all;
  
  double vm_all=0;
  double vacm = (double)vac_made;
  
  if ( nprocs > 1 ) {
	MPI_Allreduce(&naccept_double,&naccept_double_all,1,MPI_DOUBLE,MPI_SUM,world);	
	MPI_Allreduce(&vacm,&vm_all,1,MPI_DOUBLE,MPI_SUM,world);
  }
  else {
	naccept_double_all = naccept_double;
	vm_all = vacm;
  }
  if (solve) sprintf(strtmp,"%10g %12.0Lf %14d %<10d %12.0Lf",time,naccept_double_all,0,0,vm_all);
  else {
    double nattempt_double = (double)nattempt;
    double nattempt_double_all;
    if ( nprocs > 1 ) {
		MPI_Allreduce(&nattempt_double,&nattempt_double_all,1,MPI_DOUBLE,MPI_SUM,world);
	}
	else {
		nattempt_double_all = nattempt_double;
	}	
    sprintf(strtmp,"%10g %12.0Lf %14.0Lf %10d %12.0Lf",
	    time,naccept_double_all,nattempt_double_all-naccept_double_all,nsweeps,vm_all);
  }
//  check_state();
//	double dum = count_vacant();
  vac_made=0;	
}

/////////////////////////////////////////// CHECKING LATTICE /////////////////////////////////////////////////

/* ----------------------------------------------------------------------
 count pore sites
 ------------------------------------------------------------------------- */

double AppSinter::count_vacant()
{// just checking for errors
	double vacant_sites = 0; 
	
	for ( int i = 0; i < nlocal; ++i )
		if ( spin[i] == VACANT )
			vacant_sites++;
	double vacant_sites_all;
	MPI_Allreduce(&vacant_sites, &vacant_sites_all, 1, MPI_DOUBLE, MPI_SUM, world);
	
	if ( me == 0 )
		printf("time: %lf empty sites: %Lg\n", time, vacant_sites_all);
	
	return vacant_sites_all;		
}

/* ----------------------------------------------------------------------
 count grain sites
 ------------------------------------------------------------------------- */

double AppSinter::count_grain_sites()
{// just checking for errors
	double grain_sites = 0;
	
	for ( int i = 0; i < nlocal; ++i )
		if ( spin[i] != FRAME && spin[i] != VACANT )
			grain_sites++;
	double grain_sites_all;
	MPI_Allreduce(&grain_sites, &grain_sites_all, 1, MPI_DOUBLE, MPI_SUM, world);
	
//	if ( me == 0 )
//		printf("time: %lf grain sites: %d\n", time, grain_sites_all);
	
	return grain_sites_all;		
}

/* ----------------------------------------------------------------------
 calculate density
 ------------------------------------------------------------------------- */


double AppSinter::calculate_density()
{

	double grain_sites = 0;
	double total_sites = 0;
	
	int xgrid, ygrid, zgrid;

	for ( int i = 0; i < nlocal; i++ ) {
		global_to_grid( id[i], xgrid, ygrid, zgrid );
		if (	(xgrid > xstart_density) && (xgrid < xend_density) &&
				(ygrid > ystart_density) && (ygrid < yend_density) &&
				(zgrid > zstart_density) && (zgrid < zend_density) ) {
			total_sites++;
			if ( spin[i] > VACANT ) {
				grain_sites++;
			}
		}
	}
	vector<double> local_density_info( 2 );
	local_density_info[0] = grain_sites;
	local_density_info[1] = total_sites;
	
//	printf( "grain sites: %d total sites: %d\n", (int)grain_sites, (int)total_sites );
	
	vector<double> density_info_all( 2, 0 );
	
	MPI_Allreduce(&local_density_info[0], &density_info_all[0], 2, MPI_DOUBLE, MPI_SUM, world);
	
	double density = density_info_all[0] / density_info_all[1];
	
	return density;

}
