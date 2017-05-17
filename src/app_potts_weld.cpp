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

#include <stdio.h>
#include <cstdlib>
#include <cstddef>
#include <string>
#include <cstring>
#include <cmath>
#include <set>
#include <vector>
#include "error.h"
#include "domain.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>

#include "app_potts_weld.h"
#include "weld_geometry.h"
#include "teardrop.h"

using namespace SPPARKS_NS;
using weld::pool_shape::ShapeType;
using std::vector;

namespace DEMO {

   enum InchesPerMinute {I=20,II=25,III=30};

   vector< vector<double> > get_demo_control_points(InchesPerMinute ipm) {
       int n=5;
       vector< vector<double> > cp(n);
       for(int i=0;i<n;i++)
          cp[i]=vector<double>(2);
       double a[]={0.0,0.0}, b[]={0.0,0.0}, c[]={0.0,0.0}, d[]={0.0,0.0}, e[]={0.0,0.0,};
       switch(ipm){
          case I:
             a[0]=0.00;   a[1]=-6.35;
             b[0]=1.50;   b[1]=a[1];
             c[0]=3.00;   c[1]=-3.25;
             d[0]=2.75;   d[1]=0.00;
             e[0]=0.00;   e[1]=d[1];
             break;
          case II:
             a[0]=0.00;   a[1]=-7.35;
             b[0]=1.00;   b[1]=a[1];
             c[0]=2.60;   c[1]=-3.75;
             d[0]=2.80;   d[1]=0.00;
             e[0]=0.00;   e[1]=d[1];
             break;
          case III:
             a[0]=0.00;   a[1]=-7.35;
             b[0]=0.70;   b[1]=a[1];
             c[0]=2.00;   c[1]=-3.75;
             d[0]=2.475;  d[1]=0.00;
             e[0]=0.00;   e[1]=d[1];
             break;
       }

       cp[0][0]= a[0];  cp[0][1]= a[1];
       cp[1][0]= b[0];  cp[1][1]= b[1];
       cp[2][0]= c[0];  cp[2][1]= c[1];
       cp[3][0]= d[0];  cp[3][1]= d[1];
       cp[4][0]= e[0];  cp[4][1]= e[1];
       return cp;
   }

}

/* ---------------------------------------------------------------------- */

AppPottsWeld::AppPottsWeld(SPPARKS *spk, int narg, char **arg) : 
   AppPotts(spk,narg,arg), yp(0.0), alpha(0.5), beta(0.75), velocity(12500.0), 
   haz(-1.0),  distance(nullptr),
   random_park(std::atof(arg[2])), simulation_time(0.0),
   pulse_amplitude(0.0), pulse_step_frequency(1.0), shape_type(ShapeType::undefined), 
   width(-1.0), length(-1.0), teardrop_control_points()


{
   // only error check for this class, not derived classes
   if (std::strcmp(arg[0],"potts/weld") != 0 || narg != 7 )
      error->all(FLERR,"Illegal app_style in 'potts/weld' command");

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
   alpha =     std::atof(arg[3]); // relative size of pool shape at bottom compared to top
   beta =      std::atof(arg[4]); // defines curvature of pool shape through thickness 
   velocity =  std::atof(arg[5]); // travel speed: units=lattice spacings per unit time
   haz =       std::atof(arg[6]); // distance which defines 'heat affected zone'

   // add the double array
   recreate_arrays();

}

/* ---------------------------------------------------------------------- */

void AppPottsWeld::print_potts_weld_params() const {

   std::cout << "nspins: " << std::setw(7) << nspins << std::endl;
   std::cout << "pool position: " << std::setw(5) << yp << std::endl;
   std::cout << "alpha: " << std::setw(7) << alpha << std::endl;
   std::cout << "beta: " << std::setw(7) << beta << std::endl;
   std::cout << "haz: " << std::setw(7) << haz << std::endl;

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
      if (pulse_amplitude < 0.0 || pulse_amplitude >= 1.0) error->all(FLERR,"Illegal pulse command");
      pulse_step_frequency = atoi(arg[1]);
      if (pulse_step_frequency <= 2) error->all(FLERR,"Illegal pulse command");
   } else if(strcmp(command,"weld_shape_ellipse")==0) {
      if(narg != 2) error->all(FLERR,"Illegal 'weld_shape_ellipse' command");
      shape_type=ShapeType::ellipse;
      width =  std::atof(arg[0]); // weld pool width: ellipse minor axis along x-direction
      length = std::atof(arg[1]); // weld pool length: ellipse major axis along y-direction
   } else if(strcmp(command,"weld_shape_teardrop")==0) {
      shape_type=ShapeType::teardrop;
      if(strcmp(arg[0],"width")==0){ 
        width =  std::atof(arg[1]); // weld pool width (x-coordinate axis)
      } else error->all(FLERR,"Unrecognized 'weld_shape_teardrop' command.  Expected 'width'.");
      if(strcmp(arg[2],"case")==0){ 
         if(strcmp("I",arg[3])==0) teardrop_control_points=DEMO::get_demo_control_points(DEMO::InchesPerMinute::I);
         else if(strcmp("II",arg[3])==0) teardrop_control_points=DEMO::get_demo_control_points(DEMO::InchesPerMinute::II);
         else if(strcmp("III",arg[3])==0) teardrop_control_points=DEMO::get_demo_control_points(DEMO::InchesPerMinute::III);
         else if(strcmp("control_points",arg[3])==0) error->all(FLERR,"NOT IMPLEMENTED 'weld_shape_teardrop control_points' command.");
         else error->all(FLERR,"Unrecognized 'weld_shape_teardrop case' command");
      } else error->all(FLERR,"Unrecognized 'weld_shape_teardrop' command");
   } 
   else error->all(FLERR,"Unrecognized AppPottsWeld::input_app(char* command, ...) command.");
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

   // std::cout << "AppPottsWeld::user_update; dt =  " << dt << std::endl;
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

   weld::pool_shape::PoolShape *shape=nullptr;
   switch(shape_type) {
      case ShapeType::ellipse:
         shape=new weld::pool_shape::EllipticBezier(p*width,p*length,thickness,alpha,beta);
         break;
      case ShapeType::teardrop:
         shape=new weld::pool_shape::Teardrop(thickness,p*width,alpha,beta,teardrop_control_points);
         break;
   }

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
      double d=shape->distance(xo);
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

   delete shape;
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
