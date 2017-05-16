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

#ifndef SPK_AM_ELLIPSOID_H
#define SPK_AM_ELLIPSOID_H

#include <math.h>

namespace RASTER {

namespace pool_shape {

   class AmEllipsoid : public PoolShape {

      private:
         // Parameters 
         // Const is good!
         const double _p_width, _p_depth, _tail_length, _cap_height, _HAZ_width, _tail_HAZ;

      public:
         AmEllipsoid(double p_width, double p_depth, double tail_length, double cap_height, double HAZ_width, double tail_HAZ):
         	_p_width(p_width), _p_depth(p_depth), _tail_length(tail_length), _cap_height(cap_height), _HAZ_width(HAZ_width), _tail_HAZ(tail_HAZ) {}

         virtual ~AmEllipsoid(){}
         
         
         bool is_inside(const double *XYZ) const {
        	double active_p_length;
        	
        	//If z is positive, we're above the current layer and can't be inside
        	if (XYZ[2] > 0) return false;
        	
        	//Change active_p_length depending on whether XYZ[1] is pos or neg
			 if (XYZ[0] < 0) {
			 	active_p_length = _tail_length;	
			 }
			 else {
			 	active_p_length = _cap_height;
			 }		        	
            
            return pow(XYZ[1]/(_p_width * 0.5),2) + pow(XYZ[0]/active_p_length,2) + pow(XYZ[2]/_p_depth,2) < 1.0;
         }


         virtual double distance(const double *XYZ) const {

			 double active_p_length;
			 int vector_length;
			 double d = -1.0;
			 double working_tail_haz = _tail_HAZ - _tail_length;

			 			 
			 //Test to see if we're inside the melt pool, return -1 if so
			 if(is_inside(XYZ)) return d;
			 
			 //Change active_p_length depending on whether XYZ[1] is pos or neg
			 if (XYZ[0] < 0) {
			 	active_p_length = _tail_length;	
			 }
			 else {
			 	active_p_length = _cap_height;
			 }			
			 
			 //Determine the length from the melt pool center to our point.
			 //This will be our maximum possible iteration
			 vector_length = sqrt(XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1] + XYZ[2] * XYZ[2]);
			 	

			 
			//Calculate distance from trailing ellipsoid.
			//We should be able to return distance for any coordinate (although this would be expensive for far away points)
			for(int i = 0; i < vector_length; i++) {
			
				if (pow(XYZ[1]/(_p_width * 0.5 + i),2) + pow(XYZ[0]/(active_p_length + i),2) + pow(XYZ[2]/(_p_depth + i),2) <= 1) {
					
					//Smallest shell that the point is inside
					d = i;
					break;
				}
			}
							 	
            return d;
        }
        
        //This is probably the better place to do all the checks that we're currently doing in app_update
        
   };

}

}

#endif
