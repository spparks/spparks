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

#include "string.h"
#include "math.h"
#include "app_potts_weld_jom.h"
#include "random_park.h"
#include "error.h"
#include "domain.h"
#include "lattice.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsWeldJOM::AppPottsWeldJOM(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  // only error check for this class, not derived classes

  if (strcmp(arg[0],"potts/weld_jom") == 0 && narg != 10 )
    error->all(FLERR,"Illegal app_style command");
  nspins = atof(arg[1]);
  Wwidth = atoi(arg[2]);
  Wlength = atoi(arg[3]);
  Wcap = atoi(arg[4]);
  Haz = atoi(arg[5]);
  StartWeld = atoi(arg[6]);
  vel = atof(arg[7]);
  weld_type = atoi(arg[8]);
  exp_factor = atof(arg[9]);

  if(weld_type < 1 || weld_type > 5)
    error->all(FLERR,"Illegal app_style command weld_type is not correctly designated");
    
  //Lets specify some defaults for the ellipsoid parameters
  ellipsoid_depth = 0; //0.33 * domain->boxzhi/domain->lattice->zlattice;
  deep_width = 0; //0.25 * domain->boxxhi/domain->lattice->xlattice;
  deep_length = 0; //0.25 * domain->boxyhi/domain->lattice->ylattice;
    
  allow_kmc = 0;
  
  //Add the mobility array
  ndouble = 1;
  recreate_arrays();
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */
void AppPottsWeldJOM::grow_app()
{
  spin = iarray[0];
  MobilityOut = darray[0];
}

/* ----------------------------------------------------------------------
   Define additional input commands for the weld app
------------------------------------------------------------------------- */
void AppPottsWeldJOM::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"deep_width") == 0) {
  	if (narg != 1) error->all(FLERR, "Illegal deep_width command (provide one positive number)");
  	deep_width = atoi( arg[0]);
  	if(deep_width < 0) error->all(FLERR, "Illegal deep_width command (cannot be negative)");
  }
  
  else if (strcmp(command, "deep_length") == 0 ){
  	if (narg != 1) error->all(FLERR, "Illegal deep_length command (provide one positive number)");
  	deep_length = atoi( arg[0]);
  	if(deep_length < 0) error->all(FLERR, "Illegal deep_length command (cannot be negative)");
  }
  
  else if (strcmp(command, "ellipsoid_depth") == 0) {
  	if (narg != 1) error->all(FLERR, "Illegal ellipsoid_depth command (provide one positive number)");
  	ellipsoid_depth = atoi(arg[0]);
  	if(ellipsoid_depth < 0) error->all(FLERR, "Illegal ellipsoid_depth command (cannot be negative)");
  }
  
  else error->all(FLERR,"Unrecognized command");
  
}


/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsWeldJOM::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  dt_sweep = 1.0/maxneigh;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (spin[i] < 1 || spin[i] > nspins) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
  
  if(ellipsoid_depth == 0) {
  	ellipsoid_depth = 0.33 * domain->boxzhi/domain->lattice->zlattice;
  }
  if(deep_width == 0) {
  	deep_width = 0.25 * domain->boxxhi/domain->lattice->xlattice;
  }
  if(deep_length == 0) {
  	deep_length = 0.25 * domain->boxzhi/domain->lattice->zlattice;
  }
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPottsWeldJOM::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy(i);
  int y_weld;
  double Weld_width;
  double n = log10(Wwidth)/log10(Wcap);
  double n2 = log10(Wwidth)/log10(Wlength-Wcap);
  double k = pow((double)(Wlength-Wcap),0.5)/Wwidth;
  double Mobility = 1.0;
  double MobilityShallow = 0.0;
  double MobilityDeep = 0.0;

    
  // Check if weld spot is one of the traditional shapes
    if(weld_type > 2 && weld_type < 6 ) {

      if(time > StartWeld){
         Mobility = 0.0;
         y_weld = (int) (vel * (time - StartWeld));
         int y = xyz[i][1];
         if(y >= y_weld && y <= y_weld+Wlength){
             
            if(y >= y_weld + Wlength - Wcap)
             
               Weld_width = Wwidth - pow(y-y_weld-Wlength+Wcap, n);
             
            else {// if y >= y_weld and y < y_weld + Wlength
                
               if(weld_type == 3) Weld_width = Wwidth/(Wlength - Wcap)*(y-y_weld);
               else if(weld_type == 4) Weld_width = Wwidth - pow(Wlength - Wcap - y + y_weld, n2);
               else if(weld_type == 5) Weld_width = pow((y-y_weld),0.5)/k;
            }
           

        double m = (2.0/(double)(Haz - Weld_width));        

        int x = xyz[i][0];
            nx = domain->boxxhi/domain->lattice->xlattice;
        if(x > (nx-Haz)/2 && x <= nx/2)
            Mobility = m * (x - (nx-Haz)/2);
        else if(x > nx/2 && x < (nx - (nx-Haz)/2))
            Mobility = 1 - m * (x - (nx+Weld_width)/2);
        if(Mobility > 1.0 || Mobility < 0.0){ 
           Mobility = 0.0;
           spin[i] = (int) (nspins*random->uniform());
        }
         }
      }
    }
    
    
    /*We can use the shallow melt pool parameters to only do that and create a "Goldak"-esque double ellipsoid melt pool*/
    else if(weld_type == 1) {
    
    	int WorkingHAZ = (int)(Haz - Wwidth)/2.0;
        
        if (time > StartWeld) {
        
            nx = domain->boxxhi/domain->lattice->xlattice;
            ny = domain->boxyhi/domain->lattice->ylattice;
            nz = domain->boxzhi/domain->lattice->zlattice;
            
            Mobility = 0.0;
            y_weld = (int) (vel * (time - StartWeld));
            int x = xyz[i][0];
            int y = xyz[i][1];
            int z = xyz[i][2];
            
            
            // Define the wide & long, shallow melt ellipsoid
            if (y >= (y_weld - Wlength - WorkingHAZ) && y <= y_weld + deep_length/2.0 && x >= (nx/2 - Haz/2) && x <= (nx/2 + Haz/2)) {
                
                if (pow(((x - nx/2)/(Wwidth * 0.5)),2) + pow(((y - (y_weld - Wlength/2.0))/(Wlength * 0.5)),2) + pow(((z - nz)/(float)ellipsoid_depth),2) <= 1)
                {
                    Mobility = 1.1;
                }
                else {
                    
                    for (int k = 0; k <= WorkingHAZ; k++) {
                        if (pow((x - nx/2)/(Wwidth*0.5 + k),2) + pow(((y - (y_weld - Wlength/2.0))/(Wlength * 0.5 + k)),2) + pow(((z - nz)/((float)ellipsoid_depth + k)),2) <= 1 && MobilityShallow < exp(-k * exp_factor))
                        {
                            Mobility = exp(-k * exp_factor);
                        }
                    }
                    
                }
            }
			if(Mobility > 1.0 || Mobility < 0.0){
                spin[i] = (int) (nspins*random->uniform());
            }
        }
    }
    
    /* Make a keyhole shaped melt pool. The keyhole shape is made of two ellipsoids.
     One shallow and wide one, and another smaller one that goes through the height of the plate.
     It will be important to keep the parameters for each straight. */
    else if(weld_type == 2) {
        
        int WorkingHAZ = (int)(Haz - Wwidth)/2.0;
        
        if (time > StartWeld) {
        
            nx = domain->boxxhi/domain->lattice->xlattice;
            ny = domain->boxyhi/domain->lattice->ylattice;
            nz = domain->boxzhi/domain->lattice->zlattice;
            
            Mobility = 0.0;
            y_weld = (int) (vel * (time - StartWeld));
            int x = xyz[i][0];
            int y = xyz[i][1];
            int z = xyz[i][2];
            
            
            // Define the wide & long, shallow melt ellipsoid
            if (y >= (y_weld - Wlength - WorkingHAZ) && y <= y_weld + deep_length/2.0 && x >= (nx/2 - Haz/2) && x <= (nx/2 + Haz/2)) {
                
                if (pow(((x - nx/2)/(Wwidth * 0.5)),2) + pow(((y - (y_weld - Wlength/2.0))/(Wlength * 0.5)),2) + pow(((z - nz)/(float)ellipsoid_depth),2) <= 1)
                {
                    MobilityShallow = 1.1;
                }
                else {
                    
                    for (int k = 0; k <= WorkingHAZ; k++) {
                        if (pow((x - nx/2)/(Wwidth*0.5 + k),2) + pow(((y - (y_weld - Wlength/2.0))/(Wlength * 0.5 + k)),2) + pow(((z - nz)/((float)ellipsoid_depth + k)),2) <= 1 && MobilityShallow < exp(-k * exp_factor))
                        {
                            MobilityShallow = exp(-k * exp_factor);
                        }
                    }
                    
                }
            	//Set the center of the deep ellipsoid 1/4 of the shallow pool's length from the leading edge
                int y_deep = y_weld - Wlength/4.0;
                
                //Now lets define the deep ellipsoid
                if (pow(((x - nx/2)/(deep_width*0.5)),2) + pow(((y - y_deep)/(deep_length * 0.5)),2) + pow(((z - nz)/(nz* 1.1)),2) <= 1) {
                        MobilityDeep = 1.1;
                    }
            
                else {
                        
                    // We want the deep HAZ to trail the melt pool some, but we don't want the HAZ to be larger than the pool, so lets just
                    // create a gradient within the size of the regular pool
                    for (int k = 0; k <= deep_width/2.0; k++) {
                        if (pow((x - nx/2)/((deep_width*0.5 + k)),2) + pow(((y - (y_deep - Wlength/8.0))/(deep_length * 0.5 + k)),2) + pow((z - nz)/((nz* 1.1 + k)),2) <= 1 && MobilityDeep < 0.5 * exp(-k * exp_factor))
                        {
                            //We want the trailing mobility proifle to have a similar shpae to the main one, but a smaller overall value. Lets define a
                            //Coeffiecent that reduces the value
                            MobilityDeep = 0.5 * exp(-k * exp_factor);
                        }
                    }
                }
            }
        
            // Find which ellipsoid has the greatest mobility at the site
            if (z > nz - ellipsoid_depth && MobilityDeep > 1) {
                Mobility = MobilityDeep;
            }
            else {
                
                if (MobilityShallow > MobilityDeep) {
                    Mobility = MobilityShallow;
                }
                else {
                    Mobility = MobilityDeep;
                }
            }
            
            if(Mobility > 1.0 || Mobility < 0.0){
                spin[i] = (int) (nspins*random->uniform());
            }
        }
    }


  // events = spin flips to neighboring site different than self

    int j,m,value;
    int nevent = 0;
    
    
    
    if((Mobility > 0.0) && (Mobility < 1.0))    //(spin[i] != nspins) another criteria to exclude gg interaction
    {
        for (j = 0; j < numneigh[i]; j++) {
            value = spin[neighbor[i][j]];
            if (value == spin[i] || value == nspins) continue;
            for (m = 0; m < nevent; m++)
            if (value == unique[m]) break;
            if (m < nevent) continue;
            unique[nevent++] = value;
        }
        
        if (nevent == 0) return;
        int iran = (int) (nevent*random->uniform());
        if (iran >= nevent) iran = nevent-1;
        spin[i] = unique[iran];
        double efinal = site_energy(i);
        
        // accept or reject via Boltzmann criterion
        
        if (efinal <= einitial) {
            if (random->uniform() > Mobility){
                spin[i] = oldstate;
            }
        } else if (temperature == 0.0) {
            spin[i] = oldstate;
        } else if (random->uniform() > Mobility * exp((einitial-efinal)*t_inverse)) {
            spin[i] = oldstate;
        }
        
        //random is a uniform random distribution between 0 and 1
        //effectively giving all spins within the molten pool a random spin between 1 and 100
        //so nspins should be much larger than 100 (e.g. >= 10,000
        
        ////////////////////else if (Mobility > 1.0) spin[i] = (random->uniform() * 100+1);  //do i need this?
    }
    MobilityOut[i] = Mobility;
    
    
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


