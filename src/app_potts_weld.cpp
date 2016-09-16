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

/* ----------------------------------------------------------------------
   Contributing authors: John Mitchell (Sandia)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstddef>
#include <string>
#include <cstring>
#include <cmath>
#include <set>
#include "error.h"
#include "domain.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>

#include "app_potts_weld.h"
#include "weld_geometry.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsWeld::AppPottsWeld(SPPARKS *spk, int narg, char **arg) : 
   AppPotts(spk,narg,arg), yp(2500.0), width(2000.0), length(3000.0),
   alpha(0.5), beta(0.75), velocity(12500.0), haz(1000.0),  distance(0),
   random_park(std::atof(arg[2])), simulation_time(0.0),
   pulse_amplitude(0.0), pulse_step_frequency(1.0)
{
   // only error check for this class, not derived classes
   if (std::strcmp(arg[0],"potts/weld") != 0 && narg != 7 )
      error->all(FLERR,"Illegal app_style command");

   // Flag which forces 'callback' to this app each step time 'time' is updated; 
   // See 'allow_app_update' in app_lattice.h
   allow_app_update=1;
   // Adding 'distance' array; number of 'double' values per site
   // app.h
   ndouble = 1;

   // app_potts.cpp
   nspins = atoi(arg[1]);

   // app_potts_weld model parameters
   yp  =       std::atof(arg[2]); // initial pool position along y-axis
   width =     std::atof(arg[3]); // weld pool width: ellipse minor axis along x-direction
   length =    std::atof(arg[4]); // weld pool length: ellipse major axis along y-direction
   alpha =     std::atof(arg[5]); // relative size of pool shape at bottom compared to top
   beta =      std::atof(arg[6]); // defines curvature of pool shape through thickness 
   velocity =  std::atof(arg[7]); // travel speed: units=lattice spacings per unit time
   haz =       std::atof(arg[8]); // distance which defines 'heat affected zone'

   //print_potts_weld_params();

   // add the double array
   recreate_arrays();

}

/* ---------------------------------------------------------------------- */

void AppPottsWeld::print_potts_weld_params() const {

   std::cout << "nspins: " << std::setw(7) << nspins << std::endl;
   std::cout << "pool position: " << std::setw(5) << yp << std::endl;
   std::cout << "pool width: " << std::setw(7) << width << std::endl;
   std::cout << "pool length: " << std::setw(7) << length<< std::endl;
   std::cout << "alpha: " << std::setw(7) << alpha << std::endl;

}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsWeld::grow_app()
{
  spin = iarray[0];
  distance = darray[0];
}

/* ---------------------------------------------------------------------- */

void AppPottsWeld::init_app()
{

   // Print status
   // if(0==domain->me){
   //    if(screen) fprintf(screen,"Running potts/weld init_app()...\n");
   //    if(logfile) fprintf(logfile,"Running potts/weld init_app()...\n");
   // }

   delete [] sites;
   delete [] unique;
   sites = new int[1 + maxneigh];
   unique = new int[1 + maxneigh];
   dt_sweep = 1.0/maxneigh;

   int flag = 0;
   for (int i = 0; i < nlocal; i++){
      if (spin[i] < 0 || spin[i] > nspins) flag = 1;
   }
   int flagall;
   MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
   if (flagall) error->all(FLERR,"One or more sites have invalid values");

   // Initialize pool position and distances
   this->app_update(0.0);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppPottsWeld::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"pulse") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse command");
    pulse_amplitude = atof(arg[0]);
    if (pulse_amplitude < 0.0 || pulse_amplitude >= 1.0) 
      error->all(FLERR,"Illegal pulse command");
    pulse_step_frequency = atoi(arg[1]);
    if (pulse_step_frequency <= 2)
      error->all(FLERR,"Illegal pulse command");
  } else error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppPottsWeld::app_update(double dt)
{

   // Accumulate simulation time
   simulation_time+=dt;
   /*
    * START: Pulse calculations
    * NOTE: 'p' is used below to scale size of elliptical pool
    */
   double t = simulation_time;
   // Represents background power -- essentially defined by 
   //     original ellipsoidal dimensions
   double b=1.0;
   // Amplitude of pulse
   double A=pulse_amplitude;
   // Frequency of pulse
   double omega=2*M_PI/pulse_step_frequency;
   // Scalar pulse 'p' which scales size of ellipsoid
   double p=b+A*std::sin(omega*t);
   /*
    * END: Pulse calculations
    */

   //std::cout << "AppPottsWeld::user_update; dt =  " << dt << std::endl;
   // num site along each axis
   const Domain *d=this->domain;

   // domain dimension; ASSUMES lattice constant is same for each dimension
   double lx=std::abs(d->xprd);
   double lz=std::abs(d->zprd);
   double thickness=lz;

   // Update pool position
   // Location of weld pool is a function of time along y-axis but is fixed along x,z axes
   // xp: x-position of pool is at center of domain along x-axis
   // yp: y-position of pool moves along y-axis with 'velocity'; 
   // zp: z-position of pool is relative to top surface of weld pool
   double xp=lx/2.0;
   yp+=velocity*dt;
   double zp=thickness;

   weld::pool_shape::EllipticBezier eb(p*width,p*length,thickness,alpha,beta);

   // Lattice point location relative to 'pool' position
   double xo[]={0,0,0};

   // For each site i:
   //    Compute 'closest' distance from observation point to weld pool surface
   for(int i=0;i<nlocal;i++){

      // site coordinates (spparks coordinate system)
      // xyz -- public member 'app.h'

      xo[0]=xyz[i][0]-xp;
      xo[1]=xyz[i][1]-yp;
      xo[2]=xyz[i][2]-zp;

      // Compute distance for this site
      // This function returns d=-1 for points within pool;
      // Otherwise it returns the distance associated with the closest point
      double d=eb.distance(xo);
      distance[i]=d;

      /* -------------------------------------------------------------------
      For sites within pool, assign random spins so that when they exit the
      pool, they will be singleton spins and generally not part of a grain.
      ---------------------------------------------------------------------- */
      if(d<0){

         // Random number between [1,nspins] inclusive
         int ran = random_park.irandom(nspins);
         spin[i]=ran;
      }
   }
}

/* --------------------------------------------------------
  This function would make a nice 'lambda' in the future.
----------------------------------------------------------- */

double AppPottsWeld::compute_mobility(int site) const {

   double d=distance[site];
   //   mobility = 0 @ d<0
   //   mobility = 1 @ d=0 pool/solid interface
   //   mobility = 0 @ d>=haz distance defining heat effected zone
   //   Otherwise, mobility decreases linearly for 0<=d<haz
   double mobility=1-d/haz;
   if (d<0 || d>=haz)
      mobility=0.0;
   return mobility;
}

/* ----------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------- */

void AppPottsWeld::site_event_rejection(int i, RandomPark *random){

   int oldstate = spin[i];
   double einitial = site_energy(i);
   double mobility=compute_mobility(i);
   int nevent = 0;

   if(mobility>0.0){
      /*
       * Finding unique neighbor spins
       */
      std::set<int> my_unique;
      for (int j = 0; j < numneigh[i]; j++) {
         int value = spin[neighbor[i][j]];
         double j_distance = distance[neighbor[i][j]];
         // Don't flip to 'spins' on interior since they are not 'grains'
         // distance <= 0 are 'pinned' sites on the interior of weld pool;
         if (value == spin[i] || j_distance < 0.0) continue;
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
      } else if (0.0==temperature) {
         spin[i] = oldstate;
      } else if (random->uniform() > mobility * exp((einitial-efinal)*t_inverse)) {
         spin[i] = oldstate;
      }
  }

  if (spin[i] != oldstate) naccept++;
}
