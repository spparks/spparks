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
   Contributing author: Theron Rodgers and John Mitchell (Sandia)
------------------------------------------------------------------------- */

#include "string.h"
#include "math.h"
#include "random_park.h"
#include "solve.h"
#include "lattice.h"
#include "domain.h"
#include "error.h"
#include <string>
#include <iostream>
#include <limits>
#include <array>
#include <numeric>
#include "app_potts_am_bezier.h"

using std::numeric_limits;
using std::array;
using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsAmBezier::
AppPottsAmBezier(SPPARKS *spk, int narg, char **arg) :
   PottsAmPathParser(spk,narg,arg),  pool_width(0), pool_depth(0),haz(0),
   xt{0,0,0,0,0},yt{0,0,0,0,0},zs{0,0,0,0,0},
   beta_y(1.0), beta_z(0.5), distance(nullptr),
   bezier(), HAZ(), d_table()
{

   // only error check for this class, not derived classes
   if (strcmp(arg[0],"potts/am/bezier") != 0 || narg != 5 )
      error->all(FLERR,"Illegal app_style command");

   nspins        = std::stoi(arg[1]); //Number of spins
   pool_width    = std::stod(arg[2]); 
   pool_depth    = std::stod(arg[3]); 
   haz           = std::stod(arg[4]); //distance which defines 'heat affected zone'

   //Define the layer object, this might work better in init_app
   ndouble = 1;
   allow_app_update = 1;
   recreate_arrays();
}

/* ----------------------------------------------------------------------
   Define additional input commands for the AM app
------------------------------------------------------------------------- */

void AppPottsAmBezier::input_app(char *command, int narg, char **arg)
{
   /*
   parse top and spine curve control points
   potts/am/bezier control_points x x1 x2 x3 x4 x5
   potts/am/bezier control_points <y or z>  y2 y3 y4 
   */
   if (strcmp(command,"am") == 0) {
      parse_am(narg,arg);
   } else if (strcmp(command,"potts/am/bezier") == 0) {
    parse_am_bezier(narg,arg);
   } else error->all(FLERR,"Unrecognized command");
}

void AppPottsAmBezier::parse_am_bezier(int narg, char** arg)
{
   if (strcmp(arg[0],"control_points") == 0) {
      // printf("%s\n", " PottsAmPathParser::input_app \'am_path\'\n");
      if(strcmp(arg[1],"x")==0){
         if (7!=narg) error->all(FLERR,"Illegal 'potts/am/bezier' command. 'control_points x' requires 5 args.");
         for(int i=2;i<narg;i++)
            xt[i-2]=std::stod(arg[i]);
      } else if(strcmp(arg[1],"y")==0){
         if (5!=narg) error->all(FLERR,"Illegal 'potts/am/bezier' command. 'control_points y' requires 3 args.");
         yt[0]=0;
         yt[4]=0;
         for(int i=2;i<narg;i++)
            yt[i-1]=std::stod(arg[i]);
      } else if(strcmp(arg[1],"z")==0){
         if (5!=narg) error->all(FLERR,"Illegal 'potts/am/bezier' command. 'control_points z' requires 3 args.");
         zs[0]=0;
         zs[4]=0;
         for(int i=2;i<narg;i++)
            zs[i-1]=std::stod(arg[i]);
      } else {error->all(FLERR,"Illegal 'potts/am/bezier' command. Expected keyword 'x or y or z'");}
   } else if (strcmp(arg[1],"beta") == 0) {
      if (3!=narg) error->all(FLERR,"Illegal 'potts/am/bezier' command. 'beta' requires 2 args.");
      beta_y=std::stod(arg[1]);
      beta_z=std::stod(arg[2]);
   }

}
/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppPottsAmBezier::grow_app()
{
  spin = iarray[0];
  distance = darray[0];
}


/*
 * Function used by 'dLookup'
 */
array<size_t,3>
AppPottsAmBezier::get_ijk(size_t index, size_t nx, size_t ny) const {
   // computes i,j,k triplet from an index computed as follows
   //    index=i+j*nx+k*nx*ny
   // index goes fastest along x, then y, then finally z
   // Like loop over points along x, then go to anothyer

   size_t nxny=nx*ny;
   size_t k=std::floor(index/nxny);
   size_t r=index-k*nxny;
   size_t j=floor(r/nx);
   size_t i=r-j*nx;
   return std::move(array<size_t,3>{i,j,k});

}

/*
 * Function used by 'dLookup'; each proc calls this function 
 * to compute its portion of distance table.
 */
std::vector<double>
AppPottsAmBezier::
compute_dLookup_table
(int num_points, 
 size_t i0, 
 size_t j0, 
 size_t k0,
 const std::vector<double>& dx,
 const std::vector<double>& dz
)const {

   /*
    * Each proc computes distance for set of points defined by i0,j0,k0;
    * 'num_points' defines number of points for which distance is computed 
    * on each proc.  Calling function defines these inputs.
    * arg: i0,j0,k0 define starting point in table for each proc;
    *      Table is 3d;
    *      Inner iteration is on x (fastest), then y, then z (slowest);
    *      Inputs are computed by caller of this function; 
    *      i0 -- x index starting point
    *      j0 -- y index starting point
    *      k0 -- z index starting point
    * arg: num_points is number of points in table to be computed; computed 
    *      points are stored contiguously in table 
    * arg: dx, dz defines coordinates for computation of table; dy varies over [0,ny]; 
    *      table is computed for each point x,y,z formed by tensor product
    *      of dx dy dz; distance is computed for each point in table. x-components
    *      are iterated fastest, then y, then z (slowest) per loops below.
    */
   // Dimensions of d_table
   size_t nx=dim_d_table[0];
   size_t ny=dim_d_table[1];
   size_t nz=dim_d_table[2];
   std::vector<double> my_table(num_points);
   size_t count=0;
   for(size_t k=k0;k<nz;k++){
      k0=0;
      double z=dz[k];
      for(size_t j=j0;j<ny;j++){
         j0=0;
         double y=static_cast<double>(j);
         for(size_t i=i0;i<nx;i++){
            i0=0;
            double x=dx[i];
            // Position of point relative to pool
            array<double,3> r{x,y,z};
            Point point_r(r[0],r[1],r[2]);
            // If point is outside HAZ -- don't waste time 
            //    computing distance
            double d;
            if(HAZ.contains(point_r)){
               d=bezier.distance(r);
            } else {
               d=haz;
            }
            if(d>=haz)
               d=haz;
            my_table[count]=d;
            count++;
            // only compute distance for 'num_points' on each proc
            if(num_points==count) {
               return std::move(my_table);
            }
         }
      }
   }
   return std::move(my_table);
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */
void AppPottsAmBezier::init_app()
{
   // M is number of ranks
   int rank, M;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &M);

   // print control points
   //printf("Control points\n");
   //printf("\tx={%10.6f, %10.6f, %10.6f, %10.6f, %10.6f}\n",xt[0],xt[1],xt[2],xt[3],xt[4]);
   //printf("\ty={%10.6f, %10.6f, %10.6f, %10.6f, %10.6f}\n",yt[0],yt[1],yt[2],yt[3],yt[4]);
   //printf("\tz={%10.6f, %10.6f, %10.6f, %10.6f, %10.6f}\n",zs[0],zs[1],zs[2],zs[3],zs[4]);

   // Initialize bezier pool geometry
   bezier=AM_Bezier(xt,yt,zs,beta_y,beta_z,pool_width,pool_depth);
   double length=bezier.get_pool_length();
   double width=bezier.get_pool_width();
   // Print final pool dimensions
   if(0==rank){
      double depth=bezier.get_pool_depth();
      printf("void AppPottsAmBezier::init_app(); melt pool dimensions\n");
      printf("\tpool length=%7.1f\n",length);
      printf("\tpool width=%7.1f\n",width);
      printf("\tpool depth=%7.1f\n",depth);
   }

   array<double,3> c=bezier.get_pool_bb_center();

   {
      // Compute rectangular HAZ zone around pool
      Point low (c[0]-length/2-haz,c[1]-pool_width/2-haz,c[2]-pool_depth/2-haz);
      Point high(c[0]+length/2+haz,c[1]+pool_width/2+haz,c[2]+pool_depth/2+haz);
      HAZ=RASTER::rectangular_range(low,high);
   }

   // Run base class init_app
   // 'build_layer_z' variable in base class is initialized
   // Need to pass pool center z-coordinate so that its included
   init_app_am();

   /*
    * Create dLookup table for storing pre-computed distances
    * +1 on each guarantees table extent will be larger than HAZ; 
    * Only points within table are looked up; this helps to ensure
    * that table lookups do not get index out of bounds error
    */
   size_t nx=static_cast<size_t>(std::ceil(length+2*haz)+1);
   size_t ny=static_cast<size_t>(std::ceil(pool_width/2+haz)+1);
   size_t nz=static_cast<size_t>(std::ceil(pool_depth+haz)+1);
   // Make sure z-dimension of domain is sufficient 
   if(nz>(domain->boxzhi-domain->boxzlo))
      nz=domain->boxzhi-domain->boxzlo;
   // Ensure odd # points along x; will have equal number of positive 
   // and negative points with zero @ center
   // CAUTION here; doing 'double' arithmetic with integers
   double xlo=(-(double)(nx-1)/(double)2);
   double zlo=-(double)nz;
   if(0==nx%2) nx+=1;
   // +1 to include end points for z
   nz+=1;

   dim_d_table[0]=nx;
   dim_d_table[1]=ny;
   dim_d_table[2]=nz;
   std::vector<double> dx(nx,0), dz(nz,0);
   std::iota(dx.begin(),dx.end(),xlo);
   std::iota(dz.begin(),dz.end(),zlo);
   std::reverse(dz.begin(),dz.end());
   // N is number of points stored in table
   size_t N=nx*ny*nz;
   // Table stores number of points each processor will compute
   // 'f' is number of points all procs will compute
   size_t f=static_cast<size_t>(std::floor(N/M));
   // Number of points computed on each proc
   std::vector<int> np(M,f);
   // Number of points left over distributed evenly to first 'r' procs
   size_t r=N%M;
   for(size_t i=0;i<r;i++)
      np[i]+=1;

   // Check that N==sum(np[i], i=0,M-1) across all processors M
   if (std::accumulate(np.begin(),np.end(),0) != N)
      error->all(FLERR,"dLookUp error; sum(num points) across processors should equal N\n");

   // Initialize table of distances to haz
   d_table.resize(N,haz);

   // Index ptr into table used for MPI communication of table;
   // Each proc computes a set of points; rdispl stores pointer 
   // into table for each proc; 
   std::vector<int> rdispl(M,0);
   // Each proc will compute a subset of points from set of N points
   // The subset of points computed on each proc has a starting 
   // i,j,k index
   // Compute starting i,j,k indices for each point in table
   // Also assign ptr index for MPI communication below
   std::vector<size_t> start(3*M,0);
   size_t ptr=0;
   for(size_t p=0;p<M;p++){
      // Index ptr into table for each proc used for MPI communication
      rdispl[p]=ptr;
      size_t index=ptr;
      array<size_t,3> ijk=get_ijk(index,nx,ny);
      start[3*p]=ijk[0];
      start[3*p+1]=ijk[1];
      start[3*p+2]=ijk[2];
      ptr+=np[p];
   }

   // MPI all to all communication
   // Each proc will end up with a full copy of distance table
   // All to all communication forms the union of distances 
   //    calculated on each processor 
   // MPI_Alltoallv
   // https://www.open-mpi.org/doc/v4.0/man3/MPI_Alltoallv.3.php
   // https://rookiehpc.github.io/mpi/docs/mpi_alltoallv/index.html
   {
      // Each proc computes distance for its set of points
      size_t i0=start[3*rank+0];
      size_t j0=start[3*rank+1];
      size_t k0=start[3*rank+2];
      std::vector<double> sendbuff=compute_dLookup_table(np[rank],i0,j0,k0,dx,dz);
      const void *sbuff=static_cast<const void*>(sendbuff.data());
      // This processor sends same data to all processors (including itself)
      // Number of points sent is 'count'; 
      const int count=np[rank];
      std::vector<int> sendcounts(M,count);
      const int* sc=sendcounts.data();
      // Displacement/position of data for all processors is also the same 
      std::vector<int> sdispl(M,0);
      const int* sd=sdispl.data();
      // All processors have allocated full array 'd_table';
      // All processors will receive data into d_table 
      // Each processors may potentially receive different number of points -- np array
      // Each processors will have a different displacement -- rdispl array
      void *rbuff=static_cast<void*>(d_table.data());
      const int* rcounts=np.data();
      const int* rd=rdispl.data();
      MPI_Alltoallv(sbuff,sc,sd,MPI_DOUBLE,rbuff,rcounts,rd,MPI_DOUBLE,MPI_COMM_WORLD);
      {
         // Sync all processes relating to compute_dLookup_table
         int e=MPI_Barrier(MPI_COMM_WORLD);
         if (0!=e)
            error->all(FLERR,"MPI_Barrier call failed relating to sync up of dLookup setup.\n");
      }
   }

   // Compute distance function based upon initial pool position
   app_update(0.0);
}

/* ----------------------------------------------------------------------
	Use app update to set the new position of the weld pool and determine the
	mobilities for the new configuration
 ------------------------------------------------------------------------- */
void AppPottsAmBezier::app_update(double dt)
{

   // Move pool
   bool moved=app_update_am(dt);
   if(!moved){
      return;
   }

	//Use the new position as input to the distance calculation
   size_t nx=dim_d_table[0];
   size_t ny=dim_d_table[1];
   size_t nz=dim_d_table[2];

   // nx is odd; this produces an even number
   size_t cx=nx/2;

	// Go through all the local sites and lookup distance from table.
   // Bi-linear interpolate in x-y plane
   double d;
	for(int i=0;i<nlocal;i++){

		// SPPARKS lattice site
		double XYZ[]={xyz[i][0],xyz[i][1],xyz[i][2]};

		// Lattice point location relative to 'pool' position
		Point xyz_r_p=compute_position_relative_to_pool(XYZ);
      // If point is outside HAZ -- don't waste time 
      //    computing distance
      if(HAZ.contains(xyz_r_p)){
         array<double,3> xo{xyz_r_p[0],xyz_r_p[1],xyz_r_p[2]};
         // Translate point for table lookup since there are 
         //   negative x-component values for position of sites 
         //   relative to pool 
         size_t x=cx+std::floor(xo[0]);
         size_t xp1=x+1;
         // This coordinate value is used to interpolate x
         // Distance between sites is assumed to be 1.0
         double xi=xo[0]-std::floor(xo[0]);
         size_t y=std::floor(std::abs(xo[1]));
         size_t yp1=y+1;
         // This coordinate value is used to interpolate y
         // Distance between sites is assumed to be 1.0
         double eta=std::abs(xo[1])-std::floor(std::abs(xo[1]));
         size_t z=std::floor(std::abs(xo[2]));
         size_t i0=x+y*nx+z*nx*ny;
         size_t i1=xp1+y*nx+z*nx*ny;
         size_t i2=xp1+yp1*nx+z*nx*ny;
         size_t i3=x+yp1*nx+z*nx*ny;
         double d0=d_table[i0];
         double d1=d_table[i1];
         double d2=d_table[i2];
         double d3=d_table[i3];
         // bilinear interpolation
         distance[i]=d0+(d3-d0)*eta+(d1-d0+(d0-d1+d2-d3)*eta)*xi;
      } else {
         distance[i]=haz;
      }

	}

}


/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility.
 ------------------------------------------------------------------------- */

double 
AppPottsAmBezier::
compute_mobility(int site) const {

   double d=distance[site];
   //   mobility = 0 @ d<0
   //   mobility = 1 @ d=0 pool/solid interface
   //   mobility = 0 @ d>=haz distance defining heat effected zone
   //   Otherwise, mobility decreases linearly for 0<=d<haz
   double mobility=0.0;
   if (d<=0)
      mobility=1.0;
   else if(d>0 && d<haz)
      mobility=1-d/haz;
   else
      mobility=0.0;
   return mobility;

}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
void AppPottsAmBezier::site_event_rejection(int site, RandomPark *random) {
   double mobility=compute_mobility(site);

   // Inside melt pool
   // Assign random values to sites
   if (mobility <= 0)
      return;
   // Far outside melt pool beyond HAZ
   else if(mobility >= 1.0){
      spin[site] = (int) (nspins*random->uniform());
      return;
   }

   // events = spin flips to neighboring site different than self
   int oldstate = spin[site];
   double einitial = site_energy(site);
   int j,m,value;
   int nevent = 0;
   int z = xyz[site][2];

   if((mobility > 0.0) && (mobility < 1.0)) {
      for (j = 0; j < numneigh[site]; j++) {
         value = spin[neighbor[site][j]];
         double j_distance = distance[neighbor[site][j]];
         // d<0 corresponds with interior pool sites; don't flip to those
         if (value == spin[site] || j_distance<0) continue;
            for (m = 0; m < nevent; m++)
               // don't flip to a site which is same
               if (value == unique[m]) break;
            // true if didn't match existing spin value in m loop; 
            // therefore break and goto next neighbor j
            if (m < nevent) continue;
         unique[nevent++] = value;
      }

      if (0==nevent) return;
      int iran = (int) (nevent*random->uniform());
      if (iran >= nevent) iran = nevent-1;
         spin[site] = unique[iran];
      double efinal = site_energy(site);

      // accept or reject via Boltzmann criterion
      if (efinal <= einitial) {
         if (random->uniform() > mobility){
            spin[site] = oldstate;
         }
      } else if (0.0 == temperature) {
         spin[site] = oldstate;
      } else if (random->uniform() > mobility * exp((einitial-efinal)*t_inverse)) {
         spin[site] = oldstate;
      }

      if (spin[site] != oldstate) naccept++;

      // set mask if site could not have changed
      // if site changed, unset mask of sites with affected propensity
      // OK to change mask of ghost sites since never used

      if (Lmask) {
         if (einitial < 0.5*numneigh[site]) mask[site] = 1;
         if (spin[site] != oldstate)
         for (int j = 0; j < numneigh[site]; j++)
            mask[neighbor[site][j]] = 0;
      }
   }
}
