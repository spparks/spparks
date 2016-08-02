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
   Contributing authors: Eric Homer (BYU)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_phasefield.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "timer.h"

using namespace SPPARKS_NS;

// same as in create_sites.cpp and diag_cluster.cpp and lattice.cpp

enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
  FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

/* ---------------------------------------------------------------------- */

AppPottsPhaseField::AppPottsPhaseField(SPPARKS *spk, int narg, char **arg) : 
  AppPottsNeighOnly(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 1;
  allow_kmc = 0;
  allow_rejection = 1;
  allow_app_update = 1;
  allow_masking = 0;
  numrandom = 2;

  // need full neighbor lists of the 1st layer ghost sites

  delpropensity = 2;
  
  // add the double array

  recreate_arrays();

  // parse arguments for PottsPhaseField class only, not children

  if (strcmp(style,"potts/pfm") != 0) return;

  if (narg < 11) error->all(FLERR,"Illegal app_style command");
  
  // check the number of spins

  nspins = atoi(arg[1]);
  if (nspins <= 0) error->all(FLERR,"Illegal app_style command");
  if (nspins % 2)
    error->all(FLERR,"App potts/pfm must have even # of spins");
  
  phaseChangeInt = nspins/2;
  
  // dt_phasefield is set as a multiple of dt_rkmc in setup_end_app();
  // set the multiple here. ( dt_phasefield = dt_rkmc / dt_phasefield_mult )

  dt_phasefield_mult = atof(arg[2]); 
  
  // set the variables for the energetics and evolution

  gamma = atof(arg[3]);
  M_c = atof(arg[4]);
  kappaC = atof(arg[5]);
  a_1 = a_2 = atof(arg[6]);

  c_1 = atof(arg[7]);
  c_2 = atof(arg[8]);
  c_3 = atof(arg[9]);
  c_4 = atof(arg[10]);
  
  // set other default values

  nlocal_app = 0;
  cnew = NULL;
  pf_resetfield = false;
  warn_concentration_deviation = 0;
  warn_concentration_deviation_all = 0;
  cmap_ready = false;
  print_cmap = false;
  enforceConcentrationLimits = false;
  initialize_values = false;
  dimension = 0;
  cmap = NULL;

  // parse optional keywords

  int iarg = 11;

  while (iarg < narg) {
    
    if (strcmp(arg[iarg],"reset_phasefield") == 0) {
      iarg++;
      if (iarg >= narg) error->all(FLERR,"Illegal app_style potts/pfm command");
      if (strcmp(arg[iarg],"yes") == 0) pf_resetfield = true;
      else if (strcmp(arg[iarg],"no") == 0) pf_resetfield = false;
      else error->all(FLERR,"Illegal app_style potts/pfm command");
    } else if (strcmp(arg[iarg],"print_connectivity") == 0) {
      iarg++;
      if (iarg >= narg) error->all(FLERR,"Illegal app_style potts/pfm command");
      if (strcmp(arg[iarg],"yes") == 0) print_cmap = true;
      else if (strcmp(arg[iarg],"no") == 0) print_cmap = false;
      else error->all(FLERR,"Illegal app_style potts/pfm command");
    } else if (strcmp(arg[iarg],"initialize_values") == 0) {
      iarg++;
      if (iarg >= narg) error->all(FLERR,"Illegal app_style potts/pfm command");
      if (strcmp(arg[iarg],"yes") == 0) initialize_values = true;
      else if (strcmp(arg[iarg],"no") == 0) initialize_values = false;
      else error->all(FLERR,"Illegal app_style potts/pfm command");
    } else if (strcmp(arg[iarg],"enforce_concentration_limits") == 0) {
      iarg++;
      if (iarg >= narg) error->all(FLERR,"Illegal app_style potts/pfm command");
      if (strcmp(arg[iarg],"yes") == 0) {
	enforceConcentrationLimits = true;
	warn_concentration_deviation = 1;
	warn_concentration_deviation_all = 1;
      } else if (strcmp(arg[iarg],"no") == 0) {
	enforceConcentrationLimits = false; 
	warn_concentration_deviation = 0;
	warn_concentration_deviation_all = 0;
      } else error->all(FLERR,"Illegal app_style potts/pfm command");
    } else error->all(FLERR,"Illegal app_style potts/pfm command");
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

AppPottsPhaseField::~AppPottsPhaseField()
{
  if (pf_resetfield) {
    memory->sfree(pf_resetlist);
    memory->sfree(pf_resetlistvals);
  }
  if (nlocal_app) {
    memory->sfree(cnew);
    nlocal_app=0;
  }
  if (cmap)
    memory->sfree(cmap);
}

void AppPottsPhaseField::init_app() 
{
  // init the parent class first

  AppPottsNeighOnly::init_app();
  
  // setup the connectivity map

  if (!cmap_ready) setup_connectivity_map();
  
  // initialize values if command is called

  if (initialize_values) init_values();
  
  // setup array for easy reset

  if (pf_resetfield) set_phasefield_resetfield();

}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::setup_end_app() 
{
  // set the time step for phase field

  dt_phasefield = dt_rkmc / dt_phasefield_mult;
  
  // add code here at a later date to check the stability of the phase field
}

/* ----------------------------------------------------------------------
 set site value ptrs each time iarray/darray are reallocated
 ------------------------------------------------------------------------- */

void AppPottsPhaseField::grow_app()
{
  // set integer pointers

  spin = iarray[0];
  phase = iarray[1];
  
  // setup the initial phase values

  for (int i=0; i<nlocal+nghost; i++)
    set_site_phase(i);
  
  // setup double pointers

  c = darray[0];
  
  // grow cnew locally so it doesn't have to be communicated

  if (nlocal_app < nlocal) {
    nlocal_app = nlocal;

    cnew = (double *) 
    memory->srealloc(cnew,nlocal_app*sizeof(double),"app_potts_pf:cnew");
  }
}

/* ----------------------------------------------------------------------
 perform finite difference on a single  site
 ------------------------------------------------------------------------- */

void AppPottsPhaseField::site_event_finitedifference(int i)
{

  int j,jj;
  
  int ispin = spin[i];
  double qi_alpha,qi_beta,qs_alpha,qs_beta;
  
  double C_sum[3],qAlpha_sum[3],qBeta_sum[3];
  
  for (j=0; j<dimension; j++) {
    C_sum[j]=0.0;
    qAlpha_sum[j]=0.0;
    qBeta_sum[j]=0.0;
  }

  if (phase[i]==0) {
    qi_alpha=1.0;
    qi_beta=0.0;
  }
  else {
    qi_alpha=0.0;
    qi_beta=1.0;
  }
  
  // set values used for convenience

  double factor = 4.0*gamma + 2.0*a_1*qi_alpha + 2.0*a_2*qi_beta;
  double C=c[i];
  double CmC3 = C-c_3;
  double CmC4 = C-c_4;
  
  // the following code will automatically perform 1-,2-, or 3-D
  // finite difference, central in space and backward in time

  double val=0.0;
  
  // perform finite difference on all the cells in the list

  for (j=0; j<2*dimension; j++) {
    
    // site

    int s=neighbor[i][cmap[j]];

    // phase

    if (phase[s]==0) {
      qs_alpha=1.0;
      qs_beta=0.0;
    }
    else {
      qs_alpha=0.0;
      qs_beta=1.0;
    }
    
    // calculate contribution from D sum

    val += factor*c[s];
    val += 2.0*a_1*CmC3*qs_alpha + 2.0*a_2*CmC4*qs_beta;
    val += -1.0*kappaC * (-4.0*c[s] + c[neighbor[s][cmap[j]]]);
    
    // calculate the sign for the first order derivatives

    double sign;
    if (j%2) {
      sign=+1.0;
      jj=(j-1)/2;
    }
    else {
      sign=-1.0;
      jj=j/2;
    }
    
    // calculate the contribution from C, alpha and beta sum

    C_sum[jj] += sign * c[s];
    qAlpha_sum[jj] += sign * qs_alpha;
    qBeta_sum[jj] += sign * qs_beta;
    
  }
  
  // add the contribution from C, alpha and beta sums

  for (j=0; j<dimension; j++)
    val += C_sum[j] * (a_1*qAlpha_sum[j] + a_2*qBeta_sum[j]);
  
  
  // add the contribution from the central position

  val += -2.0 * dimension * 
    (4.0*gamma + 4.0*a_1*qi_alpha + 4.0*a_2*qi_beta + 3.0*kappaC)*C;
  val += 4.0 * dimension * 
    (a_1*c_3*qi_alpha + a_2*c_4*qi_beta);
    
  cnew[i] = C + (dt_phasefield * M_c) * val;
  
  if (enforceConcentrationLimits) {
    if (cnew[i] > 1.0) cnew[i]=1.0;
    if (cnew[i] < 0.0) cnew[i]=0.0;
  }
  
  if (!warn_concentration_deviation && (cnew[i] > 1.0 || cnew[i] < 0.0))
    warn_concentration_deviation=1;  
}

/* ----------------------------------------------------------------------
 compute energy of site
 ------------------------------------------------------------------------- */

double AppPottsPhaseField::site_energy(int i)
{

  // get the energy without the gradient term

  double energy = site_energy_no_gradient(i);

  // add the gradient energy for accurate energy reporting
   
   double val=0.0;

   // perform finite difference on all the cells in the list

   for (int j=0; j<2*dimension; j++) {
   
     // site
     int s=neighbor[i][cmap[j]];
     
     if (j%2==0) 
       val -= c[s];
     else
       val += c[s];
   }
   
   energy += 2.0 * kappaC * pow(val / 2.0,2.0);
  
  return energy;
}

/* ----------------------------------------------------------------------
 compute energy of site without the gradient term for efficient site event rejection
 ------------------------------------------------------------------------- */

double AppPottsPhaseField::site_energy_no_gradient(int i)
{
  int isite = spin[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != spin[neighbor[i][j]]) eng++;
  
  
  // now add the energy from deltaF 

  double energy = (double) eng;
  double C=c[i];
  
  energy+= (gamma)*( pow(C-c_1,2.0) + pow(c_2-C,2.0) );

  // now add the energy from the alpha or beta phase

  if (phase[i]==0) {
    energy+= (a_1)*pow(C-c_3,2.0);
  } else {
    energy+= (a_2)*pow(c_4-C,2.0);
  }
    
  /*--------------------------------------------------
   The gradient energy term is not calculated here for 
   computational efficiency.  This function is called 
   by site_event rejection where the gradient term 
   won't change and therefore this is more 
   computationally efficient.  The real site energy 
   term is calculated in site_energy.
  --------------------------------------------------*/
  
  return energy;
}

/* ----------------------------------------------------------------------
 rKMC method
 perform a site event with no null bin rejection
 flip to random neighbor spin without null bin
 technically this is an incorrect rejection-KMC algorithm
 ------------------------------------------------------------------------- */

void AppPottsPhaseField::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy_no_gradient(i);
  
  // events = spin flips to neighboring site different than self
  
  int j,m,value;
  int nevent = 0;
  
  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == spin[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }
  
  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform());
  if (iran >= nevent) iran = nevent-1;
  spin[i] = unique[iran];
  set_site_phase(i);
  double efinal = site_energy_no_gradient(i);
  
  // accept or reject via Boltzmann criterion
  
  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    spin[i] = oldstate;
    set_site_phase(i);
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    spin[i] = oldstate;
    set_site_phase(i);
  }
  
  if (spin[i] != oldstate) naccept++;
}

/* ----------------------------------------------------------------------
 iterate through the phase field solution 
 ------------------------------------------------------------------------- */

void AppPottsPhaseField::app_update(double stoptime)
{
  double localtime=0.0;
  
  // communicate all sites to make sure it's up-to-date when it starts

  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);
  
  int done = 0;
  while (!done) {

    // iterate through all the sets

    for (int i=0; i<nlocal; i++)
      site_event_finitedifference(i);
    
    // copy updated phase field into old concentration field

    for (int i=0; i<nlocal; i++)
      c[i]=cnew[i];
    
    // reset the certain values if appropriate

    if (pf_resetfield) {
      for (int i=0; i < pf_nresetlist; i++) 
        c[pf_resetlist[i]]=pf_resetlistvals[i];
    }
    timer->stamp(TIME_SOLVE);
    
    // re-sync all the data

    comm->all();
    
    // check for concentration devation warnings

    if (!warn_concentration_deviation_all) {
      MPI_Allreduce(&warn_concentration_deviation,
                    &warn_concentration_deviation_all,1,
                    MPI_INT,MPI_SUM,world);
      if (warn_concentration_deviation_all && me ==0) {
        if (screen) fprintf(screen,
            "Warning: Concentration has deviated outside [0,1] range\n");
        if (logfile) fprintf(logfile,
            "Warning: Concentration has deviated outside [0,1] range\n");
      }
    }
    
    timer->stamp(TIME_COMM);

    // increment time and determine when to finish

    localtime += dt_phasefield;
    if (localtime >= (stoptime- 1e-6)) done = 1;
    
    // throw exception when PF has iterated too far

    if (localtime > (stoptime + 1e-6)) {
      char errorstr[80];
      sprintf(errorstr,
              "PF step time (%f) in excess of Potts step time (%f)",
              stoptime,localtime);
      error->all(FLERR,errorstr);
    }
  } 
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::setup_connectivity_map()
{
  int i,j;
 
  // this check is redundant but I'm leaving it anyway

  if (domain->lattice->nbasis > 1)
    error->all(FLERR,
      "only single basis units are allowed for app_style potts/pfm");
  
  // set the dimension variable

  dimension = domain->dimension;
  
  if (!cmap)
    cmap = (int *) 
      memory->smalloc(2*dimension*sizeof(int),"app_potts_pf:cmap");
  
  
  // setup the connectivity map for the appropriate style

  int style=domain->lattice->style;
  
  if (style == LINE_2N) {
    cmap[0]=0; // x-1
    cmap[1]=1; // x+1
  }
  else if (style == SQ_4N) {
    cmap[0]=0; // x-1
    cmap[1]=3; // x+1
    cmap[2]=1; // y-1
    cmap[3]=2; // y+1
  }
  else if (style == SQ_8N) {
    cmap[0]=1; // x-1
    cmap[1]=6; // x+1
    cmap[2]=3; // y-1
    cmap[3]=4; // y+1
  }
           
  else if (style == SC_6N){
    cmap[0]=0; // x-1
    cmap[1]=5; // x+1
    cmap[2]=1; // y-1
    cmap[3]=4; // y+1
    cmap[4]=2; // z-1
    cmap[5]=3; // z+1
  }
  else if (style == SC_26N) {
    cmap[0]=4; // x-1
    cmap[1]=21; // x+1
    cmap[2]=10; // y-1
    cmap[3]=15; // y+1
    cmap[4]=12; // z-1
    cmap[5]=13; // z+1
  }
  else
    error->all(FLERR,
      "Lattice style not compatible with app_style potts/pfm");
    
  // connectivity map defined above should not change unless
  // create_sites changes. However, the following is an error check
  // to ensure that the connectivity map is correct.
  
  double lim[3][2];
  lim[0][0]=domain->boxxlo;
  lim[0][1]=domain->boxxhi;
  lim[1][0]=domain->boxylo;
  lim[1][1]=domain->boxyhi;
  lim[2][0]=domain->boxzlo;
  lim[2][1]=domain->boxzhi;
  
  int done,id=-1;
  for (i=0; i<nlocal; i++) {
    done=0;
    for (j=0; j<dimension; j++)
      if (xyz[i][j]!=lim[j][0] && xyz[i][j]!=lim[j][1])
        done++;
    
    if (done==dimension) {
      id=i;
      break;
    }
  } 
  
  if (id==-1) 
    error->all(FLERR,
      "Error checking connectivity map for app_style potts/pfm");
  
  // now check the ith site because it's not on any of the boundaries

  double pos;
  for (i=0; i<dimension; i++) {
    for (j=0; j<2; j++) {
      if (j==0)
        pos=xyz[id][i]-1;
      else 
        pos=xyz[id][i]+1;
      int s=neighbor[id][cmap[2*i+j]];
      if (pos != xyz[s][i]) {
        print_connectivity_map();
        error->all(FLERR,
          "Invalid connectivity map for app_style potts/pfm");
      }
    }
  }
  
  cmap_ready=true;
  
  // print the connectivity map if it has been asked for

  if (cmap_ready && print_cmap && me==0) print_connectivity_map();  
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::print_connectivity_map()
{
  char str[256];
  char *strptr;
  int i,j;

  strptr=str;

  if (dimension==1)
    sprintf(strptr,"  +/-     x\n");
  else if (dimension==2)
    sprintf(strptr,"  +/-     x     y\n");
  else
    sprintf(strptr,"  +/-     x     y     z\n");
  strptr += strlen(strptr);

  // cycle through the appropriate dimensions

  for (j=0; j<2; j++) {
    if (j==0)
      sprintf(strptr,"    -");
    else
      sprintf(strptr,"    +");
    strptr += strlen(strptr);
    
    for (i=0; i<dimension; i++) {
      sprintf(strptr," %5d",cmap[2*i+j]);
      strptr += strlen(strptr);
    }
    sprintf(strptr,"\n");
    strptr += strlen(strptr);
  }
    
  if (screen)
    fprintf(screen,
            "Connectivity map for app_style: potts/pfm\n%s\n",str);
  if (logfile)
    fprintf(logfile,
            "Connectivity map for app_style: potts/pfm\n%s\n",str);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::set_site_phase(int i)
{
  if ( (spin[i] - phaseChangeInt) > 0 )
    phase[i]=0; // Alpha
  else
    phase[i]=1; // Beta
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::set_phasefield_resetfield()
{

  // get the count
  
  pf_nresetlist=0;
  for (int i=0;i<nlocal;i++)
    if ( ( xyz[i][0] <= domain->boxxlo ) ||
	 ( xyz[i][0] >= (domain->boxxhi - 1) ) )
      pf_nresetlist++;
  
  // setup lists
  
  pf_resetlist = (int *) memory->
    smalloc(pf_nresetlist*sizeof(int),
	    "app_potts_pfm:pf_resetlist");
  pf_resetlistvals = (double *) memory->
    smalloc(pf_nresetlist*sizeof(double),
	    "app_potts_pfm:pf_resetlistvals");
  
  // set values in the list
  
  int counter=0;
  for (int i=0;i<nlocal;i++) {
    if ( xyz[i][0] <= domain->boxxlo ) {
      pf_resetlist[counter]=i;
      pf_resetlistvals[counter++]=0.0;
    } else if ( xyz[i][0] >= (domain->boxxhi - 1) ) {
      pf_resetlist[counter]=i;
      pf_resetlistvals[counter++]=1.0;
    }
  }
  if (counter != pf_nresetlist)
    error->all(FLERR,"Problem setting up pfsolver reset field lists");
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseField::init_values()
{
  double q_alpha,q_beta,cAlpha,cBeta;
  
  // phase=0 - alpha (a_1,c_3)

  q_alpha = 1.0; q_beta = 0.0;
  cAlpha=(gamma * (c_1 + c_2) + a_1*c_3*q_alpha + a_2*c_4*q_beta )
    / (2.0*gamma + a_1*q_alpha + a_2*q_beta);
  
  // phase=1 - beta  (a_2,c_4)

  q_alpha = 0.0; q_beta = 1.0;
  cBeta =(gamma * (c_1 + c_2) + a_1*c_3*q_alpha + a_2*c_4*q_beta )
    / (2.0*gamma + a_1*q_alpha + a_2*q_beta);
  
  for (int i=0;i<nlocal;i++) {
    set_site_phase(i);
    if (phase[i]==0)
      c[i] = cAlpha;
    else
      c[i] = cBeta;
  }

  // ensure this doesn't get called again

  initialize_values=false;
}
