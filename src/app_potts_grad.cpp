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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "app_potts_grad.h"
#include "random_park.h"
#include "domain.h"
#include "stdlib.h"
#include "lattice.h"
#include "error.h"

#include <set>
#include <limits>

using namespace SPPARKS_NS;

//enums consistent with lattice.cpp
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

/* ---------------------------------------------------------------------- */

AppPottsGrad::AppPottsGrad(SPPARKS *spk, int narg, char **arg) : 
   AppPotts(spk,narg,arg), T0(0.0), 
   grad_x(0.0), grad_y(0.0), grad_z(0.0),
   max_T(std::numeric_limits<double>::min()),
   max_M(0.0), T(0), m0(1.0), M(0), activation_energy(1.0), convert(1.0)
{
  // parse arguments for PottsGran class only, not children

  if (strcmp(style,"potts/grad") != 0) return;

   if (narg < 9 || narg > 10)
     error->all(FLERR,"Invalid amount of AppPottsGrad args");

   int n = strlen(arg[1]) + 1;
   gradient_choice = new char[n];
   strcpy(gradient_choice,arg[1]);

   nspins = atof(arg[2]);

   m0 = atof(arg[3]); // as in m = m0 * exp(-q/kT) for temp grad and initial mobility at the center for mobility grad
   if(m0 <= 0)
	   error->all(FLERR,"Invalid mobility argument");

   convert = atof(arg[4]);
   if(convert <= 0)
	   error->all(FLERR,"Invalid convert argument");

   activation_energy = atof(arg[5]);
   if(activation_energy < 0)
	   error->all(FLERR,"Invalid activation energy argument");

   T0 = atof(arg[6]);
   if(T0 < 0)
	   error->all(FLERR,"Invalid temperature argument");

   grad_x = atof(arg[7]) * convert;

   grad_y = atof(arg[8]) * convert;

   if(narg == 10){
	 grad_z = atof(arg[9]) * convert;
   }

   if(!(strcmp(gradient_choice,"temp") == 0 || strcmp(gradient_choice,"mob") == 0))
	   error->all(FLERR,"Invalid grad_style command");

   // adding temperature array 'T'
   // add mobility array 'M'
   ndouble = 2;

   //KMC not allowed for this app
   allow_kmc = 0;

   // add the double array
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsGrad::grow_app()
{
  spin = iarray[0];
  T = darray[0];
  M = darray[1];
}

void AppPottsGrad::init_app()
{
   delete [] sites;
   delete [] unique;
   sites = new int[1 + maxneigh];
   unique = new int[1 + maxneigh];

   dt_sweep = 1.0/maxneigh;


  //error checking for specific app
   int latticeStyle = domain->lattice->style;
   if(!(latticeStyle == SC_6N || latticeStyle == SC_26N || latticeStyle == SQ_4N || latticeStyle == SQ_8N)){
	   error->all(FLERR,"illegal lattice style with app");
   }


   /*
    * Initialize site temperatures;
    * Find maximum temperature across all processors;
    */
   if(strcmp(gradient_choice,"temp") == 0)
   {
	   double global_max;
	   for(int i=0;i<nlocal;i++){
		  T[i]=compute_site_temperature(i);
		  //temperature cannot be negative
		  if(T[i] < 0)
			  error->one(FLERR,"gradient caused negative temperature");

		  if(max_T<T[i]) max_T=T[i];
	   }
	   MPI_Allreduce(&max_T,&global_max,1,MPI_DOUBLE,MPI_MAX,world);
	   max_T=global_max;
	   for(int i=0;i<nlocal;i++){
		  // skip this site if pinned
		  if(M[i]<0.0) continue;
		  M[i]=compute_mobility(T[i],max_T);
	   }
   }

   /*
    * Initialize site mobilities directly from mobility gradient    *
    */
   else
   {
	   double global_mob_max;
	   for(int i=0;i<nlocal;i++){
		  // skip this site if pinned
		  if(M[i]<0.0) continue;
		  M[i]=mobility_grad(i);
		  if(max_M < M[i])
			  max_M = M[i];

		  //mobility cannot be negative
		  if(M[i] < 0)
			  error->one(FLERR,"gradient caused negative mobility");
	   }

	   MPI_Allreduce(&max_M,&global_mob_max,1,MPI_DOUBLE,MPI_MAX,world);
	   max_M = global_mob_max;

	   //scale the mobility from 0 to 1
	   for(int i = 0; i<nlocal;i++)
		   M[i]=M[i]/max_M;

   }

   int flag = 0;
   for (int i = 0; i < nlocal; i++)
      if (spin[i] < 0 || spin[i] > nspins) flag = 1;
   int flagall;
   MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
   if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

double AppPottsGrad::mobility_grad(int i){
	  double c=this->convert;
	   const Domain *d=this->domain;

	   // num site along each axis
	   int nx = domain->nx;
	   int ny = domain->ny;
	   int nz = domain->nz;

	   // domain dimension; ASSUMES lattice constant is same for each dimension
	   double lx=d->xprd;
	   double h=lx/nx;

	   // site coordinates (spparks coordinate system)
	   // xyz -- public member 'app.h'
	   double xi = xyz[i][0];
	   double yi = xyz[i][1];
	   double zi = xyz[i][2];

	   double x = xi + (1.0-nx) * h / 2.0;
	   double y = yi + (1.0-ny) * h / 2.0;
	   double z = zi + (1.0-nz) * h / 2.0;

	   return m0 + grad_x * x + grad_y * y + grad_z * z;
}

double AppPottsGrad::compute_mobility(double Ti,double max_T) const {
   
   // Ti: site temperature
   // max_T: maximum temperature across entire domain
   // Boltzman constant k = 1.3806488 × 10^{-23} Units --> (m^2 - kg / (s^2 K)) = energy/(per degree K)
   // Boltzmann constant using 'electron volts' as unit of energy
   const double K=.000086171;
   double kT=K*Ti; // k*Ti, where K is Boltzmann's constant 8.6171e-5 eV/K
   // activation energy
   double q= activation_energy;
   double max_M=m0*exp(-q/(K*max_T));

   /** 
    * The gb mobility have to be scaled to our simulations as M[i]/max_M,
    *  otherwise the simulation will run too long for no good reason.
    */
   return m0*exp(-q/(kT))/max_M;
}

double AppPottsGrad::compute_temperature(double x, double y, double z) const {
   return  T0 + grad_x * x + grad_y * y + grad_z * z;
}

double AppPottsGrad::compute_site_temperature(int i) const {
   double c=this->convert;
   const Domain *d=this->domain;

   // num site along each axis
   int nx = domain->nx;
   int ny = domain->ny;
   int nz = domain->nz;

   // domain dimension; ASSUMES lattice constant is same for each dimension
   double lx=d->xprd;
   double h=lx/nx;

   // site coordinates (spparks coordinate system)
   // xyz -- public member 'app.h'
   double xi = xyz[i][0];
   double yi = xyz[i][1];
   double zi = xyz[i][2];

   double x = xi + (1.0-nx) * h / 2.0;
   double y = yi + (1.0-ny) * h / 2.0;
   double z = zi + (1.0-nz) * h / 2.0;

   return compute_temperature(x,y,z);
}
/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsGrad::site_event_rejection(int i, RandomPark *random){

   int oldstate = spin[i];
   double einitial = site_energy(i);
   double Ti=T[i];
   double mobility=M[i];

   int nevent = 0;
   // mobilities <= 0 are 'pinned' sites
   if(mobility > 0.0){
      /*
       * Finding unique neighbor spins
       */

	 std::set<int> my_unique;
	  for (int j = 0; j < numneigh[i]; j++) {
		 int value = spin[neighbor[i][j]];
		 double j_mobility= M[neighbor[i][j]];
		 // mobilities <= 0 are 'pinned' sites
		 if (value == spin[i] || j_mobility < 0.0) continue;
		 my_unique.insert(value);
	  }
	  // how many unique spins
	  nevent=my_unique.size();

	  // scale random number according number of unique neighbor spins
	  if (nevent == 0) return;
	  int iran = (int) (nevent*random->uniform());
	  if (iran >= nevent) iran = nevent-1;


	  {
		 /*
		  * Store unique spins
		  */
		 int c=0;
		  for(std::set<int>::iterator i=my_unique.begin();i!=my_unique.end();i++){
			unique[c++]=*i;
		  }
	  }

      spin[i] = unique[iran];
      double efinal = site_energy(i);

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > mobility){
               spin[i] = oldstate;
         }
      } else if (temperature == 0.0) {
         spin[i] = oldstate;
      } else if (random->uniform() > mobility * exp((einitial-efinal)*t_inverse)) {
         spin[i] = oldstate;
      }
  }

  if (spin[i] != oldstate) naccept++;

}

AppPottsGrad::~AppPottsGrad()
{
  delete [] gradient_choice;
}
