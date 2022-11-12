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
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#include "app_potts_ca.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "timer.h"
#include "output.h"
#include "lattice.h"
#include "app.h"
#include "domain.h"

#include <vector>
#include <set>

using namespace SPPARKS_NS;

enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};

/* ---------------------------------------------------------------------- */

AppPottsCA::AppPottsCA(SPPARKS *spk, int narg, char **arg) :
  AppPottsNeighOnly(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 1;
  delpropensity = 1;
  delevent = 0;
  allow_kmc = 1;
  allow_rejection = 1;
  allow_masking = 1;
  allow_app_update = 1;
  numrandom = 1;
  
  recreate_arrays();
	
  // parse arguments
  
  if (narg != 6) error->all(FLERR,"Illegal app_style command");
	
  counter = 0;
  nspins = atoi(arg[1]);
  dt_sweep = 1.0/nspins;
	
  rc_repeats = atof(arg[2]);
  gg_repeats = atof(arg[3]);
  nuclei_pct = atof(arg[4]);
  threshold = atof(arg[5]);

  sites = unique = NULL;

  // RNG for lattice_updates()

  ranlatt = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  ranlatt->reset(seed,me,100);
}

/* ---------------------------------------------------------------------- */

AppPottsCA::~AppPottsCA()
{
  delete ranlatt;
}

/* ---------------------------------------------------------------------- 
   set site value ptrs for spins (iarray) and stored energy (darray) 
   each time iarray and darray are reallocated
------------------------------------------------------------------------- */

void AppPottsCA::grow_app()	
{ 
  spin = iarray[0];
  stored_e = darray[0];
  nuc_switch = iarray[1];
}

/* ----------------------------------------------------------------------
   Lattice Updates - following rejection_Kmc algorithms, first update each 
   site's dvalue to enact an incremental rise in stored energy and secondly 
   introduce nuclei at random locations throughout array based on a probability 
   associated with current stored energies
 
   note: lattice updating scales with number of neighbors as "nlocal" associates 
   each site with those of its neighborhood, therefore, incremental
   increases are scaled down by a prefactor: 1D(1/2), 2D(1/8), 3D(1/26).
 
   note: since some sites do not have the maximum neighbor of neighbors the
   scaling prefactor is not exact and variation in increment w.r.t. stored
   energy may be introduced at a singular time step(s)
 ------------------------------------------------------------------------- */

void AppPottsCA::app_update(double dt_caller)
{
  // determine current maximum stored energy in system (at last timestep)
  
  double current_proc_pre_eMax = 0.0;
  for (int i = 0; i < nlocal; i++)  {
    if (stored_e[i] > current_proc_pre_eMax) 
      {current_proc_pre_eMax = stored_e[i];}
  }
  double current_global_pre_eMax;
  MPI_Allreduce(&current_proc_pre_eMax, &current_global_pre_eMax, 1,
                MPI_DOUBLE, MPI_MAX, world);

  // determine dimension dependent argument(s) (e.g. aspect)

  int nz, ny, nx;
	
  if ((nz != 1) && (ny != 1) && (nx != 1)) { aspect = 3;}
	
  if ((nz == 1) && (ny != 1) && (nx != 1)) { aspect = 2;}
  if ((nz != 1) && (ny == 1) && (nx != 1)) { aspect = 2;}
  if ((nz != 1) && (ny != 1) && (nx == 1)) { aspect = 2;}
  
  if ((nz != 1) && (ny == 1) && (nx == 1)) { aspect = 1;}  
  if ((nz == 1) && (ny != 1) && (nx == 1)) { aspect = 1;} 
  if ((nz == 1) && (ny == 1) && (nx != 1)) { aspect = 1;}
  
  if ((nz == 1) && (ny == 1) && (nx == 1)) { aspect = 0;}
  
  // increase stored energies
	
  for (int i = 0; i < nlocal; i++) {
    if (stored_e[i] < 5) {
      int xcord = xyz[i][0];
      int ycord = xyz[i][1];
      int zcord = xyz[i][2];
      chance = ranlatt->irandom(MAXSMALLINT) % 100000;
      if (chance < 1000) {	
	double coeff1 = 0.0000000001 * time;
	double coeff2 = 0.0000000001 * time;
	 
        // there's a parallelization problem here
        // trying to determine stored energy based on neighbor
        // but stored energy isn't being incremented for sites between processor regions
        
	std::vector<int>::iterator it;
        std::vector<int> neighlist;
        neighlist.resize(numneigh[i]);
        for (int j = 0; j < numneigh[i]; j++)
          neighlist[j] = spin[neighbor[i][j]];
	
        std::set<int> s(neighlist.begin(), neighlist.end());
        neighlist.assign(s.begin(), s.end());

	if (stored_e[i] < 5) {	
          if (aspect == 3) stored_e[i] += ((ranlatt->irandom(MAXSMALLINT) % 10 + 1 * 0.1)/26);
          if (aspect == 2) stored_e[i] += ((ranlatt->irandom(MAXSMALLINT) % 10 + 1 * 0.1)/8);
          if (aspect == 1) stored_e[i] += ((ranlatt->irandom(MAXSMALLINT) % 1000 * 0.001)/2);
        }
      }
    }
    
    // determine current maximum stored energy in system (after energy increase)
	
    double current_proc_post_eMax = 0.0;
    for (int i = 0; i < nlocal; i++)	 {
      if (stored_e[i] > current_proc_post_eMax) {
        current_proc_post_eMax = stored_e[i];
      }
    }
    double current_global_post_eMax;
    MPI_Allreduce(&current_proc_post_eMax, &current_global_post_eMax, 1,
                  MPI_DOUBLE, MPI_MAX, world);

    // seed nuclei 
    // changed for unique spins by D.Rule
    
    for (int i = 0; i < nlocal; ++i) {
      double value = ranlatt->irandom(MAXSMALLINT) % 1000 * 0.001;
      double current_e = (stored_e[i]/current_global_post_eMax);
      
      // for ordering of nuclei
	
      if (current_e > threshold) {
        ++counter;
        if (value < nuclei_pct) {	
          if (nuc_switch[i] == 0) {
            stored_e[i] = 0;
            spin[i]= nspins + counter;
          }
        }
      }	
    }
  }
}
 
/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsCA::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  dt_sweep = 1.0/maxneigh;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (spin[i] < 1 /*|| spin[i] > nspins*/) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- 
   compute the total energy at each lattice site by either considering 
   contribution from unlike neighbors (site_energy = eng) or the locally 
   stored energy at each site (site_energy = s_eng) 
   - For grain growth (GG), site energy is totally a function of 
   surface energy from unlike neighbors (eng)
   - For nucleation events (NU), site energy is totally a function of
   stored energy at each lattice site (s_eng)
   - The argument i is a packed argument.  It is actually 2*i+C where
   C is 1 for an nucleation attempt and 0 for a gg attempt.
 ------------------------------------------------------------------------- */

double AppPottsCA::site_energy(int i)	
{ 
  // establish a int-i-based conditional operation for the selection of a 
  // site energy derived from the stored energy of the site (s_eng) or the 
  // unlike neighbors of the site (eng).

  int site = (i/2);
  int NucAttmpt = i%2;
  double total_eng = 0;
	
  // contribution due to stored energy via darray
  
  if (NucAttmpt == 1) {
    return stored_e[site];
  }
  
  // for clustered continuous nucleation the site energy could be based
  // on the number spins bordering
  // contribution due to surface energy caused by unlike neighbors
  
  // grain growth is attempted
  
  int isite = spin[site];
  int eng = 0;
  for (int j = 0; j < numneigh[site]; j++)
    if (isite != spin[neighbor[site][j]]) eng++;
  
  return (double) eng;
}

/* ----------------------------------------------------------------------
   rKMC method - Recrystallization (RC) & Grain Growth (GG)
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm . . .
   ... For (RC), recrystallized sites are updated to a stored energy of zero
   ... For (GG), flipped sites' dvalue are updated to stored energy of joined neighbor
 ------------------------------------------------------------------------- */

void AppPottsCA::site_event_rejection(int i, RandomPark *random)
{
  // notes:	
  // rc_repeats are the number of recrystallization steps
  // gg_repeats are the number of grain growth steps
  // repeat values determine a fractional proportion of RC and GG events
  //   over the entire domain imposed by comparison with a
  //   random value between [0,1]
  // defining both rc and gg values independently then determining 
  //   the ratio allow either value to be zero if desired
	
  // Recrystallization (RC) algorithm in which probability is a func (dvalue)
  // probability of RC event based on percent RC of total repeats
  
  if (1*random->uniform() <= (rc_repeats/(rc_repeats + gg_repeats))) {
    double stored_eMax = (0.00093*pow(stoptime,2)) - (0.0001*stoptime) + 1.0;

    int oldstate = spin[i];
    double einitial = site_energy(2*i+1); 
	
    // odd argument signifies that site_energy is for recrystallization step
    // tells site_energy to call the stored_e
    
    double unique_e[1 + maxneigh];

    //creates an array of 1 more than maxneigh. 2D the size is 9 doubles.
    // events = spin flips to neighboring site different than self
	
    int j,m,value;
    double s_energy;
    int nevent = 0;

    for (j = 0; j < numneigh[i]; j++) {
      if (stored_e[i] == 0) continue;

      value = spin[neighbor[i][j]];
      s_energy = stored_e[neighbor[i][j]];  
      if (stored_e[i] <= s_energy) continue; 
      for (m = 0; m <= nevent; m++)
        if (s_energy == unique_e[m]) break;
      if (m < nevent) continue;
      unique[nevent] = value;
      unique_e[nevent++] = s_energy;  
    }

    int iran = (int) (maxneigh*random->uniform());
    if (iran >= nevent) return;
	
    // neighonly spin flip/rejection	
    // update flipped site with both spin & dvalue of joined neighbor
    // says that you can update site to joined neighbors if your neighbor is a "nuclei"
    
    spin[i] = unique[iran];
    stored_e[i]	= unique_e[iran];       // or equal to zero?
		
    double efinal = stored_e[i];
	
    // accept or reject via a modified stored energy criterion

    if (efinal <= einitial) {
    } else if (temperature == 0.0) {
      spin[i] = oldstate;
      stored_e[i] = einitial;
    } else if (random->uniform() > ((fabs(efinal-einitial))/stored_eMax) ) {
      spin[i] = oldstate;
      stored_e[i] = einitial;
    }
	
    if (spin[i] != oldstate) naccept++;
	
    // set mask if site could not have changed
    // if site changed, unset mask of sites with affected propensity
    // OK to change mask of ghost sites since never used
    
    if (Lmask) {
      if (einitial < 0.5*numneigh[i]) mask[i] = 1;
      if (spin[i] != oldstate)
        for (int j = 0; j < numneigh[i]; j++)
          mask[neighbor[i][j]] = 0;
    }
  
  // Grain Growth (GG) algorithm in which probability is a func (unlike neighbors)

  } else {
    int oldstate = spin[i];
    double einitial = site_energy(2*i); // even argument signifies that site_energy is for gg step.
    double unique_e[1 + maxneigh];
	
    // events = spin flips to neighboring site different than self
	
    int j,m,value;
    double s_energy;
    int nevent = 0;
    for (j = 0; j < numneigh[i]; j++) {
      if ((spin[neighbor[i][j]]) > (nspins)) { 
        value = spin[neighbor[i][j]];
        s_energy = stored_e[neighbor[i][j]];  
        if (value == spin[i]) continue;
        for (m = 0; m < nevent; m++)
          if (value == unique[m]) break;
        if (m < nevent) continue;
        unique[nevent] = value;
        unique_e[nevent++] = s_energy;  
      }
    }
    
    if (nevent == 0) return;
    int iran = (int) (nevent*random->uniform());
    if (iran >= nevent) iran = nevent-1;
	
    // update flipped site with both spin & dvalue (dvalue = zero) of joined neighbor
    
    spin[i]		= unique[iran];
    stored_e[i]	= unique_e[iran]; 	
	
    double efinal = site_energy(2*i);
	
    // accept or reject via Boltzmann criterion
	
    if (efinal <= einitial) {
    } else if (temperature == 0.0) {
      spin[i] = oldstate;
      stored_e[i] = einitial;
    } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
      spin[i] = oldstate;
      stored_e[i] = einitial;
    }
	
    if (spin[i] != oldstate) naccept++;
	
    // set mask if site could not have changed
    // if site changed, unset mask of sites with affected propensity
    // OK to change mask of ghost sites since never used
	
    if (Lmask) {
      if (einitial < 0.5*numneigh[i]) mask[i] = 1;
      if (spin[i] != oldstate)
        for (int j = 0; j < numneigh[i]; j++)
          mask[neighbor[i][j]] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   rKMC method - Lattice Scanning 
------------------------------------------------------------------------- */

/*
void AppPottsCA::iterate_rejection(double)
{
  int i,icolor,nselect,nrange,jset;
  int *site2i;
	
  // set loop is over:
  // sectors if there are sectors and no colors
  // colors if there are colors and no sectors
  // first nsector sets if there are both sectors and colors
	
  int nset_loop = nset;
  if (bothflag) nset_loop = nsector;
	
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
					 			
      timer->stamp(TIME_SOLVE);			// signal end of serial sweep
			
      if (nprocs > 1) {
        if (sectorflag) comm->reverse_sector(iset);
        else comm->all_reverse();
        timer->stamp(TIME_COMM);		// signal end of parallel sweep
      }
    }

    // call to lattice_update()
    // increase stored energies throughout array, inttroduce nuclei into array
    
    app_update(dt_rkmc);
		
    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }	
}	
*/
