/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
   Class AppPhaseFieldThinFilm - added by Fadi Abdeljawad, fabdelj@sandia.gov
   Most of this was copied from AppPotts, AppPottsNeighOnly and PottsPhaseField.
 
   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
 
   See the README file in the top-level SPPARKS directory.
   ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_phasefield_thinfilm.h"
#include "solve.h"
#include "random_park.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

// same as in create_sites.cpp and diag_cluster.cpp and lattice.cpp

enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
     FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

/* ---------------------------------------------------------------------- */

AppPhaseFieldThinFilm::AppPhaseFieldThinFilm(SPPARKS *spk, int narg, char **arg) :
  AppPottsNeighOnly(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 104;
  allow_kmc = 0;
  allow_rejection = 1;
  allow_app_update = 1;
  allow_masking = 0;
  numrandom = 2;
  delpropensity = 2;  //need full neighbor lists of the 1st layer ghost sites
    
  // add the double array

  recreate_arrays();
    
  if (narg != 14)
    error->all(FLERR,"Illegal app_style command - PhaseField/ThinFilm");

  // set phase field variables here

  num_OP = atof(arg[2]);
  pf_dt = atof(arg[3]);
  Wphi = atof(arg[4]);
  mH = atof(arg[5]);
  kappa = atof(arg[6]);
  epsilon = atof(arg[7]);
  Lphi = atof(arg[8]);
  M_o = atof(arg[9]);
  M_vol = atof(arg[10]);
  M_eva = atof(arg[11]);
  M_sur = atof(arg[12]);
  M_gb = atof(arg[13]);
    
  // set other default values

  nlocal_app=0;
  cmap_ready=false;
  print_cmap=false;
  initialize_values=false;
  dimension=0;
  cmap=NULL;
}

/* ---------------------------------------------------------------------- */

AppPhaseFieldThinFilm::~AppPhaseFieldThinFilm()
{
  if (nlocal_app) {
    memory->destroy(chm_phi);
    nlocal_app=0;
  }
  if (cmap) memory->sfree(cmap);
}

/* ---------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"PFthinfilm/command") == 0) {
    if (narg != 2) error->all(FLERR,"Invalid command for app_style");
    else if (strcmp(arg[0],"print_connectivity") == 0) {
      if (strcmp(arg[1],"yes") == 0) print_cmap=true;
      else print_cmap=false;
    }
    else if (strcmp(arg[0],"initialize_values") == 0) {
      if (strcmp(arg[1],"yes") == 0) initialize_values=true;
      else initialize_values=false;
    }
    else
      error->all(FLERR,"Invalid command for app_style");
  } else error->all(FLERR,"Invalid command for app_style");
}

/* ---------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::init_app()
{
  //init the parent class first

  AppPottsNeighOnly::init_app();
    
  //setup the connectivity map

  if (!cmap_ready) setup_connectivity_map();
    
  //initialize values if command is called

  if (initialize_values) init_values();
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   ------------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::grow_app()
{
  // you have to keep this to trick SPPARKS into running

  spin = iarray[0]; 
    
  // setup double pointers

  c = darray[num_OP];
  GB = darray[num_OP + 1];
  chmPot = darray[num_OP + 2];
  mobility = darray[num_OP + 3];
    
  if (nlocal_app < nlocal) {
    nlocal_app = nlocal;
    memory->create(chm_phi,num_OP,nlocal_app,"app_pf_thinfilm:chm_phi");
  }
}

/* ----------------------------------------------------------------------
   perform finite difference on a single site
------------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::site_event_finitedifference1(int i)
{
  int j,jj, OP;
  int OP_start = 0;
    
  //set values used for convenience

  double C=c[i];
    
  double val2=0.0;
  double valphi1 = 0.0;
  double valphi2 = 0.0;
  double sumSquare = 0.0;
  double sumCube = 0.0;
  double four_thirds = 1.3333333333333;
    
  // calculate Sum phi^2 and Sum phi^3

  for (OP = 0; OP<num_OP; OP++) {
    double phi_val = darray[OP_start + OP][i];
    sumSquare += phi_val * phi_val;
    sumCube += phi_val * phi_val * phi_val;
  }
    
  // Finding the chemical potential of the mass density field
  // first: get the local quantities
    
  val2 = 4.0 * mH * C * (C * C - 1.0) + 
    four_thirds * Wphi * (0.5 * (C + 1.0) - 3.0 * sumSquare + 2.0 * sumCube);
    
  // perform finite difference to calculate: nabla^2 C

  for (j=0; j<2*dimension; j++) {
    // site neighbors
    int s=neighbor[i][cmap[j]];
    //calculate contribution from neighbors
    val2 += - kappa * kappa * c[s];
  }
    
  //add the contribution from the local point

  val2 += 2.0 * dimension * kappa * kappa * C;
    
  // The chemical potential of C

  chmPot[i] = val2;

  // atomic mobility of c
  // constant mobility everywhere
  //    mobility[i] = M_o;

  // surface diffusion only

  double g1_eta = 0.03125 * (C * C * C * (6.0 * C * C - 20.0) + 30.0 * C + 16.0);
  mobility[i] = M_o + M_vol * g1_eta + M_eva * (1.0 - g1_eta) + 
    M_sur * (1.0 - C) * (1.0 + C) + M_gb * 2.0 * GB[i];

  // to make sure no negative diffusivity

  if (mobility[i] < 0.0) {mobility[i] = 0.0;}
    
  // Chemical potential of each of the phi's
    
  for (OP = 0; OP<num_OP; OP++) {
    double phi_val = darray[OP_start + OP][i];
    double chmphi_val = 0.0;
        
    chmphi_val += 8.0 * Wphi * phi_val * 
      ((1 - C) - (3.0 - C) * phi_val + 2.0 * sumSquare);

      for (j=0; j<2*dimension; j++) {
        int s=neighbor[i][cmap[j]];
        chmphi_val += - epsilon * epsilon * darray[OP_start + OP][s];
      }
        
      chmphi_val += 2.0 * dimension * epsilon * epsilon * phi_val;
      chm_phi[OP][i] = chmphi_val;
    }
}

/* ----------------------------------------------------------------------
   Update rule for the C PDE and ones for phi's
------------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::site_event_finitedifference2(int i)
{
  int j,jj, OP;
    
  double dummy1 = 0.0;
  int OP_start = 0;
  double GB_field = 0.0;
    
  for (OP = 0; OP<num_OP; OP++) {
    darray[OP_start + OP][i] -= Lphi * pf_dt * chm_phi[OP][i];
    GB_field += darray[OP_start + OP][i] * (1.0 - darray[OP_start + OP][i]);
  }
  
  GB[i] = GB_field;
    
  double mob_term = 0.0;
    
  for (j=0; j<2*dimension; j++) {
    int s=neighbor[i][cmap[j]];
    mob_term += mobility[s] * chmPot[s] + mobility[i] * chmPot[s] - 
      chmPot[i] * mobility[s];
  }
    
  mob_term += -2.0 * dimension * mobility[i] * chmPot[i];
  c[i] = c[i] + 0.5 * pf_dt * mob_term;
}


/* ----------------------------------------------------------------------
   compute energy of site
   ------------------------------------------------------------------------- */

double AppPhaseFieldThinFilm::site_energy(int i)
{
  int j,jj, OP;
  int OP_start = 0;
  double four_thirds = 1.3333333333333;
    
  //set values used for convenience

  double C=c[i];
    
  double energy = mH * (C * C - 1.0) * (C * C - 1.0) + 0.33333333 * 
    Wphi * (C + 1.0) * (C + 1.0);
  double sumSquare = 0.0;
  double sumCube = 0.0;

  for (OP = 0; OP < num_OP; OP++) {
    double phi_val = darray[OP_start + OP][i];
        
    sumSquare += phi_val * phi_val;
    sumCube += phi_val * phi_val * phi_val;
  }
    
  energy += four_thirds * Wphi * (3.0 * (1.0 - C) * sumSquare - 
                                  2.0 * (3.0 - C) * sumCube + 
                                  3.0 * sumSquare * sumSquare);
  return energy;
}

/* ----------------------------------------------------------------------
   ------------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::site_event_rejection(int i, RandomPark *random)
{
  /*
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
  */
}

/* ----------------------------------------------------------------------
   iterate through the phase field solution
   ------------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::app_update(double stoptime)
{
  double localtime = 0.0;
    
  // communicate all sites to make sure up-to-date when it starts

  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);
    
  // iterate through all the sets

  for (int i=0; i<nlocal; i++) {site_event_finitedifference1(i);}
    
  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);
    
  for (int i=0; i<nlocal; i++) {site_event_finitedifference2(i);}
    
  timer->stamp(TIME_SOLVE);
    
  // re-sync all the data

  comm->all();
    
  timer->stamp(TIME_COMM);
}

/* ---------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::setup_connectivity_map()
{
  int i,j;
    
  //this check is redundant but I'm leaving it anyway

  if (domain->lattice->nbasis > 1)
    error->all(FLERR,"Only single basis units allowed for app_style PF/ThinFilm");
    
  //set the dimension variable

  dimension = domain->dimension;
    
  if (!cmap)
    cmap = (int *)
      memory->smalloc(2*dimension*sizeof(int),"app_pf_thinfilm:cmap");
    
    
  //setup the connectivity map for the appropriate style

  int style=domain->lattice->style;
    
  if (style == LINE_2N) {
    cmap[0]=0;// x-1
    cmap[1]=1;// x+1
  }
  else if (style == SQ_4N) {
    cmap[0]=0;// x-1
    cmap[1]=3;// x+1
    cmap[2]=1;// y-1
    cmap[3]=2;// y+1
  }
  else if (style == SQ_8N) {
    cmap[0]=1;// x-1
    cmap[1]=6;// x+1
    cmap[2]=3;// y-1
    cmap[3]=4;// y+1
  }
    
  else if (style == SC_6N){
    cmap[0]=0;// x-1
    cmap[1]=5;// x+1
    cmap[2]=1;// y-1
    cmap[3]=4;// y+1
    cmap[4]=2;// z-1
    cmap[5]=3;// z+1
  }
  else if (style == SC_26N) {
    cmap[0]=4; // x-1
    cmap[1]=21;// x+1
    cmap[2]=10;// y-1
    cmap[3]=15;// y+1
    cmap[4]=12;// z-1
    cmap[5]=13;// z+1
  }
  else
    error->all(FLERR,"Lattice style not compatible with app_style PF/ThinFilm");
    
  //The connectivity map defined above should not change unless
  //create_sites changes. However, the following is an error check
  //to ensure that the connectivity map is correct.
    
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
    error->all(FLERR,"Error in connectivity map for app_style PF/ThinFilm");
    
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
        error->all(FLERR,"Invalid connectivity map for app_style PF/ThinFilm");
      }
    }
  }
    
  cmap_ready=true;
    
  // print the connectivity map if it has been asked for

  if (cmap_ready && print_cmap && me==0) print_connectivity_map();
}

/* ---------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::print_connectivity_map()
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
    
  //cycle through the appropriate dimensions

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
            "Connectivity map for app_style: PF/ThinFilm\n%s\n",str);
  if (logfile)
    fprintf(logfile,
            "Connectivity map for app_style: PF/ThinFilm\n%s\n",str);
}

/* ---------------------------------------------------------------------- */

void AppPhaseFieldThinFilm::init_values()
{
  srand(685479);
    
  int OP_start = 0;
  for (int i = 0; i < nlocal; i++) {
    c[i] = 1.0;
    for (int OP = 0; OP < num_OP; OP++) 
      darray[OP_start + OP][i] = 0.5 + 0.15 * 2.0 * 
        ((double)rand() / (double)RAND_MAX - 0.5);
  }

  initialize_values=false;//ensure this doesn't get called again
}
