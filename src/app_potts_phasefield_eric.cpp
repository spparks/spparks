/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
   Class AppPottsPhaseField - added by Eric Homer, ehomer@sandia.gov
   Mar 31, 2011 - Most recent version.  Most of this was copied from 
   AppPotts and AppPottsNeighOnly.

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "app_potts_phasefield_eric.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "timer.h"
#include "table.h"

#include "style_region.h"

using namespace SPPARKS_NS;

// same as in create_sites.cpp and diag_cluster.cpp and lattice.cpp
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
  FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

#define BIG 1000;

/* ---------------------------------------------------------------------- */

AppPottsPhaseFieldEric::
AppPottsPhaseFieldEric(SPPARKS *spk, int narg, char **arg) : 
  AppPottsNeighOnly(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 2;
  allow_kmc = 0;
  allow_rejection = 1;
  allow_masking = 0;
  allow_app_update = 1;
  numrandom = 3;
  delpropensity = 2;//need full neighbor lists of the 1st layer ghost sites
  
  //add the double array

  recreate_arrays();

  nspins = atoi(arg[1]);
  if (nspins <= 0) error->all(FLERR,"Illegal app_style command");

  int nphasestemp = atoi(arg[2]);
  if (nphasestemp < 1) error->all(FLERR,"Illegal app_style command");
  
  if ((9+2*nphasestemp) != narg)
    error->all(FLERR,"Illegal app_style command - AppPottsPhaseFieldEric");
  
  //dt_phasefield is set as a multiple of dt_rkmc in setup_end_app();
  //set the multiple here. ( dt_phasefield = dt_rkmc / dt_phasefield_mult )

  dt_phasefield_mult = atof(arg[3]); 
  
  //set the variables for the energetics and evolution

  mobility=atof(arg[4]);
  ch_energy=atof(arg[5]);
  thermal_conductivity = atof(arg[6]);
  heat_of_transport = atof(arg[7]);
  specific_heat = atof(arg[8]);
  
  phases = (Phase **) 
    memory->smalloc(nphasestemp*sizeof(Phase *),"app_potts_pf:phases");
  
  int iarg=9;
  int currentspin=0;
  nphases = 0;

  for (int i=0; i<nphasestemp; i++) {
    
    int phaseID = find_phase(arg[iarg]);
    if (phaseID >= 0 ) error->all(FLERR,"Phase already exists");
    //create the new phase
    phases[nphases] = new Phase;
    
    int n = strlen(arg[iarg]) + 1;
    phases[nphases]->id = new char[n];
    strcpy(phases[nphases]->id,arg[iarg]);
    phases[nphases]->phase = nphases;
    phases[nphases]->spin_range = atoi(arg[iarg+1]);
    phases[nphases]->spin_start = currentspin+1;
    currentspin += phases[nphases]->spin_range;
    phases[nphases]->spin_end = currentspin;
    phases[nphases]->table = new Table(spk,phases[nphases]->id);
    //set other default values
    phases[nphases]->nRow = phases[nphases]->nCol = 0;
    phases[nphases]->rowvals = phases[nphases]->colvals = NULL;
    phases[nphases]->G = phases[nphases]->dGdC = NULL;
    phases[nphases]->scaleval = 1.0;
    
    iarg += 2;
    nphases++;
  }

  if (nphases != nphasestemp) error->all(FLERR,"Illegal app_style command");
  if (currentspin != nspins) error->all(FLERR,"Illegal app_style command "
                                        "Phase spins don't add to nspins");
  
  //set other default values
  nlocal_app=0;
  cnew=NULL;
  Tnew=NULL;
  pf_resetfield=false;
  pf_nresetlist = pf_maxresetlist = 0;
  warn_concentration_deviation=0;
  warn_concentration_deviation_all=0;
  cmap_ready=false;
  print_cmap=false;
  init_site_phase=false;
  enforceConcentrationLimits=false;
  dimension=0;
  latconst=0.0;
  cmap=NULL;
  
  run_nucleation = false;
  nucleation_rate = 0.0;
  
  scale_potts_energy = 1.0;
  scale_free_energy = 1.0;
  scale_comp_energy_penalty = 1.0;
  
}

/* ---------------------------------------------------------------------- */

AppPottsPhaseFieldEric::~AppPottsPhaseFieldEric()
{
  if (pf_resetfield) {
    memory->destroy(pf_resetlist);
    memory->destroy(pf_resetlistvals);
  }
  if (nlocal_app) {
    memory->sfree(cnew);
    memory->sfree(Tnew);
    nlocal_app=0;
  }
  if (cmap)
    memory->sfree(cmap);
  
  for (int i = 0; i < nphases; i++) {
    delete [] phases[i]->id;
    delete phases[i]->table;
    memory->destroy(phases[i]->dGdC);
    delete phases[i];
  }
  memory->sfree(phases);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::input_app(char *command, int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Invalid command for app_style");
  
  if (strcmp(command,"pottspf/command") == 0) {
    
    if (strcmp(arg[0],"reset_phasefield") == 0)
      set_pf_resetlist(narg,arg);
    
    else if (strcmp(arg[0],"print_connectivity") == 0) {
      if (strcmp(arg[1],"yes") == 0) print_cmap=true;
      else print_cmap=false;
    }
    
    else if (strcmp(arg[0],"set_site_phase") == 0) {
      if (strcmp(arg[1],"yes") == 0) init_site_phase=true;
      else print_cmap=false;
    }
    
    else if (strcmp(arg[0],"enforce_concentration_limits") == 0) {
      if (strcmp(arg[1],"yes") == 0) {
        enforceConcentrationLimits=true;
        //since I'm forcing the concentration limits
        //set the warn flags to 1 so that it doesn't check
        warn_concentration_deviation=1;
        warn_concentration_deviation_all=1;
        
        if (domain->me == 0) {
          if (screen) fprintf(screen,
                              "Concentration limits will be enforced\n");
          if (logfile) fprintf(logfile,
                               "Concentration limits will be enforced\n");
        }
      } 
      else {
        enforceConcentrationLimits=false; 
        //if someone turns this off reset the flags to check
        // for concentration deviations
        warn_concentration_deviation=0;
        warn_concentration_deviation_all=0;
        
        if (domain->me == 0) {
          if (screen) fprintf(screen,
                              "Concentration limits will not be enforced\n");
          if (logfile) fprintf(logfile,
                               "Concentration limits will not be enforced\n");
        }
      }
    }
    
    else if (strcmp(arg[0],"nucleation") == 0) {
      error->all(FLERR,"Nucleation command not currently available");
      run_nucleation = true;
      nucleation_rate = atof(arg[1]);
      
      if (domain->me == 0) {
        if (screen) fprintf(screen,
                            "Setting a nucleation rate of %f\n",
                            nucleation_rate);
        if (logfile) fprintf(logfile,
                             "Setting a nucleation rate of %f\n",
                             nucleation_rate);
      }
    }
    
    else if (strcmp(arg[0],"scale_energy") == 0) {
      int i=1;
      while (i < narg) {
        if (strcmp(arg[i],"potts")==0) {
          scale_potts_energy = atof(arg[i+1]);
          if (domain->me == 0) {
            if (screen) fprintf(screen,
                                "Potts interface energy will be scaled by %f\n",
                                scale_potts_energy);
            if (logfile) fprintf(logfile,
                                 "Potts interface energy will be scaled by %f\n",
                                 scale_potts_energy);
          }
        }
        else if (strcmp(arg[i],"free")==0) {
          scale_free_energy = atof(arg[i+1]);
          if (domain->me == 0) {
            if (screen) fprintf(screen,
                                "Free energy will be scaled by %f\n",
                                scale_free_energy);
            if (logfile) fprintf(logfile,
                                 "Free energy will be scaled by %f\n",
                                 scale_free_energy);
          }
        }
        else if (strcmp(arg[i],"penalty")==0) {
          scale_comp_energy_penalty = atof(arg[i+1]);
          if (domain->me == 0) {
            if (screen) fprintf(screen,
                                "Composition energy penalty will be scaled by %f\n",
                                scale_comp_energy_penalty);
            if (logfile) fprintf(logfile,
                                 "Composition energy penalty will be scaled by %f\n",
                                 scale_comp_energy_penalty);
          }
        }
        else
          error->all(FLERR,"Illegal scale_energy command");
        i+=2;
      }
    }
    else error->all(FLERR,"Invalid command for app_style");
  }
  else error->all(FLERR,"Invalid command for app_style");
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::init_app() 
{
  int i;
  //check to make sure that simulation is fully periodic
  for (i=0; i<dimension; i++)
    if (domain->periodicity[i] != 1)
      error->all(FLERR,"app_style does not support non-periodic boundaries");
  
  //init the parent class first
  AppPottsNeighOnly::init_app();
  
  //setup the connectivity map
  if (!cmap_ready) setup_connectivity_map();
  
  //process phases
  process_phases();
  
  if (init_site_phase)
    for (i=0; i<nlocal+nghost; i++)
      set_site_phase(i);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::process_phases()
{
  int i,j,k;
  
  //check to make sure all phase tables have been imported
  for (i=0; i<nphases; i++) {
    Phase *p = phases[i];
    Table *t = p->table;
    if (!t->getTableReady())
      error->all(FLERR,"Not all tables loaded for phases");
    
    //setup pointers to table values
    t->getTableDims(&p->nRow,&p->nCol);
    p->G = t->getTablePtr();
    p->rowvals = t->getRowPtr();
    p->colvals = t->getColPtr();
    if (p->dGdC == NULL)
      memory->create(p->dGdC,p->nRow,p->nCol,"app_potts_pf:p->dGdC");
    
    p->dC = p->colvals[1]-p->colvals[0];
    p->dT = p->rowvals[1]-p->rowvals[0];
    
    //colvals must be the concentration and increase from 0 to 1
    if (p->colvals[0] !=0 || p->colvals[p->nCol-1] != 1)
      error->all(FLERR,"Illegal table for phase, concentration must "
                 "increase from 0 to 1");
    for (j=1; j<p->nCol; j++)
      if (p->colvals[j] < 0 || p->colvals[j] > 1 || 
          p->colvals[j] <= p->colvals[j-1] ||
          fabs(p->colvals[j]-p->colvals[j-1] - p->dC) > 1e-6)
        error->all(FLERR,"Illegal table for phase, concentration must "
                   "increase from 0 to 1 and have even intervals");
    
    //rowvals must be the temperature and increase, error if anything less than or equal zero
    if (p->rowvals[0] <= 0)
      error->all(FLERR,"Illegal table for phase, temperature must be greater than zero");
    for (j=1; j<p->nRow; j++) {
      if (p->rowvals[j] <= 0)
        error->all(FLERR,"Illegal table for phase, temperature must be greater than zero");
      if (p->rowvals[j] <= p->rowvals[j-1] ||
          fabs(p->rowvals[j]-p->rowvals[j-1] - p->dT) > 1e-6)
        error->all(FLERR,"Illegal table for phase, temperature must increase "
                   " and have even intervals");
    }
    
    //scale the energy values if appropriate
    if (p->scaleval != scale_free_energy) {
      for (j=0; j<p->nRow; j++)
        for (k=0; k<p->nCol; k++)
          p->G[j][k] *= scale_free_energy;
      p->scaleval = scale_free_energy;
    }
    
    //calculate dGdC for the G matrix;
    double h=p->dC;
    int n=p->nCol - 1;//index of the last value
    double **G=p->G;
    double **dGdC=p->dGdC;
    for (j=0; j<p->nRow; j++) {
      //forward difference on first value
      dGdC[j][0] = (-G[j][2] + 4*G[j][1] - 3*G[j][0])/(2*h);
      //central difference on middle values
      for (int k=1; k<n; k++)
        dGdC[j][k] = (G[j][k+1] - G[j][k-1])/(2*h);
      //backward difference on last value
      dGdC[j][n] = (3*G[j][n] - 4*G[j][n-1] + G[j][n-2])/(2*h);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::setup_end_app() 
{
  //set the time step for phase field
  dt_phasefield = dt_rkmc / dt_phasefield_mult;
}

/* ----------------------------------------------------------------------
 set site value ptrs each time iarray/darray are reallocated
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::grow_app()
{
  //set integer pointers
  spin = iarray[0];
  phase = iarray[1];
  
  //setup the initial phase values
  for (int i=0; i<nlocal+nghost; i++)
    set_site_phase(i);
  
  //setup double pointers
  c = darray[0];
  T = darray[1];
  
  //grow cnew locally so it doesn't have to be communicated
  if (nlocal_app < nlocal) {
    nlocal_app = nlocal;

    cnew = (double *) 
    memory->srealloc(cnew,nlocal_app*sizeof(double),"app_potts_pf:cnew");
    Tnew = (double *) 
    memory->srealloc(Tnew,nlocal_app*sizeof(double),"app_potts_pf:Tnew");
  }
}

/* ----------------------------------------------------------------------
 perform finite difference on a single  site
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::site_event_finitedifference(int i)
{

  int p,j,jj;
  double qsum[3],fsum[3],Tsum[3],Tinvsum[3];
  
  //set initial values
  double hsquared = latconst * latconst;
  
  // the following code will automatically perform 1-,2-, or 3-D
  //finite difference, central in space and backward in time
  double valA=0.0;
  double valB=0.0;
  double valC=0.0;
  double qi, qj, fpi, fpj;
  
  for (p=0; p<nphases; p++) {
    if (phase[i] == p)
      qi = 1;
    else
      qi = 0;
    
    fpi = bilinear_interp(i,1);
    
    for (j=0; j<dimension; j++)
      qsum[j] = fsum[j] = Tsum[j] = Tinvsum[j] = 0.0;
    
    //perform finite difference on all the cells in the list
    for (j=0; j<2*dimension; j++) {
      
      //site
      int s=neighbor[i][cmap[j]];
      
      if (phase[s] == p)
        qj = 1;
      else
        qj = 0;
      
      //calculate bilinear interp at site s
      double fpj = bilinear_interp(s,1);
      
      valA += qj * fpi + qi * fpj;
      
      //calculate the sign and position for the first order derivatives
      double sign;
      if (j%2) {
        sign=+1.0;
        jj=(j-1)/2;
      }
      else {
        sign=-1.0;
        jj=j/2;
      }
      qsum[jj] += sign * qj;
      fsum[jj] += sign * fpj;
    }
    
    //add the contribution from the first order gradients
    for (j=0; j<dimension; j++)
      valA += 2 * qsum[j] * fsum[j];
    
    //add contribution from site i values
    valA -= 4.0 * dimension * qi * fpi;
  }
  
  //perform finite difference on all the cells in the list
  //for phase independent values
  for (j=0; j<2*dimension; j++) {
    
    //site
    int s=neighbor[i][cmap[j]];
    
    //calculate contribution from Cahn-Hilliard term
    valA -= (ch_energy/hsquared) * (-4.0*c[s] + c[neighbor[s][cmap[j]]]);
    
    //calculate contribution from temperature gradients
    valB += T[s];
    
    //calculate the sign and position for the first order derivatives
    double sign;
    if (j%2) {
      sign=+1.0;
      jj=(j-1)/2;
    }
    else {
      sign=-1.0;
      jj=j/2;
    }
    Tsum[jj] += sign * T[s];
    Tinvsum[jj] += sign / T[s];
    
    
  }
  //add contribution from the central position
  valA -= 6.0 * dimension * ch_energy * c[i] / hsquared;
  
  //divide by final hsquared to get the value of the Laplacian of the 
  //chemical potential
  valA /= hsquared;
  
  //add the contribution form the central position
  valB -= 2.0 * dimension * T[i];
  
  //divide by final hsquared to get the value of the Laplacian of the 
  //temperature
  valB /= hsquared;
  
  //add the contribution from the first order gradients
  for (j=0; j<dimension; j++)
    valC += Tsum[j] * Tinvsum[j];
  
  //divide by final hsquared and add the contribution from the 
  //Laplacian of the temperature
  valC = (valC / hsquared) + (valB / T[i]);
  
  
  //calculate new value of the composition
  cnew[i] = c[i] + (dt_phasefield * mobility) * 
                   (2.0 * valA + heat_of_transport * valC);
  
  //calculate new value of the temperature
  Tnew[i] = T[i] + (dt_phasefield / specific_heat) * 
                   ( (2.0 * mobility * heat_of_transport * valA ) + 
                     (thermal_conductivity * valB) );
  
  //enforce concentration limits if appropriate
  if (enforceConcentrationLimits) {
    if (cnew[i] > 1.0) cnew[i]=1.0;
    if (cnew[i] < 0.0) cnew[i]=0.0;
  }
  
  //check warnings
  if (!warn_concentration_deviation && (cnew[i] > 1.0 || cnew[i] < 0.0))
    warn_concentration_deviation=1;  
  
  if (!warn_temperature_deviation && Tnew[i] < 0.0) 
    warn_temperature_deviation=1; 
  
  //do not let Temperature go below 0
  if (Tnew[i] < 0.0) Tnew[i] = 0.0;
}

/* ----------------------------------------------------------------------
 compute energy of site
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldEric::site_energy(int i)
{
  //return the full energy value 
  return (site_energy_bulk(i) 
        + site_energy_potts(i) 
        + site_energy_phasefield(i));
}

/* ----------------------------------------------------------------------
 compute energy of site without the gradient term for efficient site event rejection
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldEric::site_energy_without_phasefield(int i)
{    
  /*------------------------------------------------------
   The phase field energy term is not calculated here for 
   computational efficiency.  This function is called by
   site_event rejection where the gradient term won't 
   change and therefore this is more computationally 
   efficient.  The real site energy term is calculated 
   in site_energy.
   ------------------------------------------------------*/
  return site_energy_bulk(i) + site_energy_potts(i);
}

/* ----------------------------------------------------------------------
 compute the potts interface energy of site
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldEric::site_energy_potts(int i)
{
  //calculate the Potts interface energy
  int isite = spin[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != spin[neighbor[i][j]]) eng++;
  
  // Now scale the energy if appropriate 
  double energy = (double) eng;
  
  return energy * scale_potts_energy;
}

/* ----------------------------------------------------------------------
 compute the bulk free energy of site
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldEric::site_energy_bulk(int i)
{
  //no sum over phases because the site only has one phase
  //energy scaling has already been accounted for
  return bilinear_interp(i,0);
}
  
/* ----------------------------------------------------------------------
 compute the bulk free energy of site
 ------------------------------------------------------------------------- */
  
double AppPottsPhaseFieldEric::site_energy_phasefield(int i)
{
  
  double val=0.0;
  //perform finite difference on all the cells in the list
  for (int j=0; j<2*dimension; j++) {
    
    int s=neighbor[i][cmap[j]];
    
    if (j%2==0) 
      val -= c[s];
    else
      val += c[s];
  }
  
  return 0.5 * ch_energy * pow(val / 2.0,2.0);
}

/* ----------------------------------------------------------------------
 rKMC method
 perform a site event with no null bin rejection
 flip to random neighbor spin without null bin
 technically this is an incorrect rejection-KMC algorithm
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy_without_phasefield(i);
  
  // events = spin flips to neighboring site different than self
  
  int j,m,value;
  int nevent = 0;

  //nucleation probability - check if running - otherwise normal site event
  if (run_nucleation && (random->uniform() < nucleation_rate))
  {

    //Generate a random spin n, irandom(n), of the opposite phase
    //For alpha nucleation sites, add phaseChangeInt

    //spin[i] = random->irandom(phaseChangeInt);
    //if (phase[i] == ALPHA) {
      //was alpha (phase=0, spin <= phaseChangeInt)
      //switch to beta (phase=1, spin > phaseChangeInt)
    //  spin[i] += phaseChangeInt;
    //}
    //else it was beta (phase=1, spin > phaseChangeInt)
    //will switch to alpha (phase=0, spin <= phaseChangeInt)
    //set_site_phase(i);
  }
  else
  {
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
    double efinal = site_energy_without_phasefield(i);
    
    // accept or reject via Boltzmann criterion
    
    if (efinal <= einitial) {
    } else if (temperature == 0.0) {
      spin[i] = oldstate;
      set_site_phase(i);
    } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
      spin[i] = oldstate;
      set_site_phase(i);
    }
  }
  
  if (spin[i] != oldstate) naccept++;
}

/* ----------------------------------------------------------------------
 iterate through the phase field solution 
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::user_update(double stoptime)
{
  double localtime=0.0;
  
  //communicate all sites to make sure it's up-to-date when it starts
  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);
  
  int done = 0;
  while (!done) {

    //iterate through all the sets
    for (int i=0; i<nlocal; i++)
      site_event_finitedifference(i);
    
    //copy updated phase field into old concentration field
    for (int i=0; i<nlocal; i++) {
      c[i]=cnew[i];
      T[i]=Tnew[i];
    }
    
    //reset the certain values if appropriate
    if (pf_resetfield) {
      for (int i=0; i < pf_nresetlist; i++) 
        c[pf_resetlist[i]]=pf_resetlistvals[i];
    }
    timer->stamp(TIME_SOLVE);
    
    //re-sync all the data
    comm->all();
    
    //check for concentration devation warnings
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
    
    //check for temperature devation warnings
    if (!warn_temperature_deviation_all) {
      MPI_Allreduce(&warn_temperature_deviation,
                    &warn_temperature_deviation_all,1,
                    MPI_INT,MPI_SUM,world);
      if (warn_temperature_deviation_all && me ==0) {
        if (screen) fprintf(screen,
                            "Warning: temperature attempted to go below 0\n");
        if (logfile) fprintf(logfile,
                             "Warning: temperature attempted to go below 0\n");
      }
    }
    
    timer->stamp(TIME_COMM);

    //increment time and determine when to finish
    localtime += dt_phasefield;
    if (localtime >= (stoptime- 1e-6)) done = 1;
    
    //throw exception when PF has iterated too far
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

void AppPottsPhaseFieldEric::setup_connectivity_map()
{
  int i,j;
 
  //this check is redundant but I'm leaving it anyway
  if (domain->lattice->nbasis > 1)
    error->all(FLERR,
      "only single basis units are allowed for app_style potts/phasefield");
  
  //set the dimension variable
  dimension = domain->dimension;
  //set the lattice constant variable
  latconst = domain->lattice->latconst;
  
  if (!cmap)
    cmap = (int *) 
      memory->smalloc(2*dimension*sizeof(int),"app_potts_pf:cmap");
  
  
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
    error->all(FLERR,
      "Lattice style not compatible with app_style potts/phasefield");
    
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
    error->all(FLERR,
      "Error checking connectivity map for app_style potts/phasefield");
  
  //now check the ith site because it's not on any of the boundaries
  double pos;
  for (i=0; i<dimension; i++) {
    for (j=0; j<2; j++) {
      if (j==0)
        pos=xyz[id][i]-latconst;
      else 
        pos=xyz[id][i]+latconst;
      int s=neighbor[id][cmap[2*i+j]];
      if (pos != xyz[s][i]) {
        print_connectivity_map();
        error->all(FLERR,
          "Invalid connectivity map for app_style potts/phasefield");
      }
    }
  }
  
  cmap_ready=true;
  
  //print the connectivity map if it has been asked for
  if (cmap_ready && print_cmap && me==0) print_connectivity_map();  
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::print_connectivity_map()
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
            "Connectivity map for app_style: potts/phasefield\n%s\n",str);
  if (logfile)
    fprintf(logfile,
            "Connectivity map for app_style: potts/phasefield\n%s\n",str);
}
  
/* ---------------------------------------------------------------------- */
  
double AppPottsPhaseFieldEric::bilinear_interp(int i,int flag)
{
  //bilinear interpolation formula from equation 13.54 in 
  //MACHINE VISION by Ramesh Jain, Rangachar Kasturi, Brian G. Schunck
  //Published by McGraw-Hill, Inc., ISBN 0-07-032018-7, 1995
  //http://www.cse.usf.edu/~r1k/MachineVisionBook/MachineVision.htm
  
  double temp = T[i];
  double comp = c[i];
  Phase *p=phases[phase[i]];
  double val;
  double **mat;
  if (flag==0) mat = p->G;
  else mat = p->dGdC;
  int iC,iT;
  double deltaC,deltaT;
  
  //The out of bounds factor is a quadratic function that is employed for 
  //  compositions outside the range of [0,1] so it will bring them back.
  //  It is a quadratic function that is valued such that additional energy at 
  //  10% above 1 or below 0 will have the value of scale_free_energy added onto
  //  the energy at 0 or 1. Essentially:
  //  G(C) = G(1 or 0) + outofboundsfactor* Cexcess^2 and
  //  dGdC(C) = 2 * outofboundsfactor * Cexcess
  double outofboundsfactor = scale_comp_energy_penalty * scale_free_energy;
  
  int outofbounds=0;
  if (comp < 0) {
    if (flag) return outofboundsfactor*2.0*comp;
    outofbounds = -1;
    iC=0;
  } else if (comp > 1) {
    if (flag) return outofboundsfactor*2.0*(comp-1.0);
    outofbounds = 1;
    iC=p->nCol-1;
  }
  if (outofbounds && flag == 0) {
    iT=binarySearch(p->rowvals,0,p->nRow-1,temp);
    if (p->rowvals[iT] > temp)
      iT--;
    if (iT == p->nRow-1)
      iT--;
    deltaT = (temp - p->rowvals[iT])/p->dT;
    val = mat[iT][iC] + deltaT*(mat[iT+1][iC] - mat[iT][iC]);
  
    if (outofbounds == -1)
      val += outofboundsfactor * comp * comp;
    else
      val += outofboundsfactor * (comp-1) * (comp-1);
    return val;
    
  }
  
  //find starting point for the matrix
  iC=binarySearch(p->colvals,0,p->nCol-1,comp);
  if (p->colvals[iC] > comp)
    iC--;
  if (iC == p->nCol-1)
    iC--;
  deltaC = (comp - p->colvals[iC])/p->dC;
  
  iT=binarySearch(p->rowvals,0,p->nRow-1,temp);
  if (p->rowvals[iT] > temp)
    iT--;
  if (iT == p->nRow-1)
    iT--;
  deltaT = (temp - p->rowvals[iT])/p->dT;
  
  val = mat[iT][iC] + deltaT*(mat[iT+1][iC] - mat[iT][iC]) + 
                      deltaC*(mat[iT][iC+1] - mat[iT][iC]) +
        deltaC*deltaT * (mat[iT][iC] - mat[iT+1][iC] -
                         mat[iT][iC+1] + mat[iT+1][iC+1]);
  
  return val;
}

/* ---------------------------------------------------------------------- */

int AppPottsPhaseFieldEric::find_phase(char *name)
{
  for (int iphase = 0; iphase < nphases; iphase++)
    if (strcmp(name,phases[iphase]->id) == 0) return iphase;
  return -1;
}
/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::set_site_phase(int i)
{
  for (int j=0; j<nphases; j++)
    if (spin[i] >= phases[j]->spin_start &&
        spin[i] <= phases[j]->spin_end)
      phase[i]=j;
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldEric::set_pf_resetlist(int narg, char **arg)
{
  if (narg < 3)
    error->all(FLERR,"Illegal pottspf/command reset_phasefield");
  
  //process the arguments one by one
  int iarg = 1;
  while (iarg < narg) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal pottspf/command "
                                  "reset_phasefield");
    
    int iregion = domain->find_region(arg[iarg]);
    if (iregion < 0) error->all(FLERR,"reset_phasefield command region ID "
                                "does not exist");
    double val = atof(arg[iarg+1]);
    
    for (int i=0;i<nlocal;i++)
      if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
        if (pf_nresetlist >= pf_maxresetlist) {
          //grow the arrays
          pf_maxresetlist += BIG;
          memory->grow(pf_resetlist,pf_maxresetlist,
                       "app_potts_pf:pf_resetlist");
          memory->grow(pf_resetlistvals,pf_maxresetlist,
                       "app_potts_pf:pf_resetlistvals");
        }
        
        //add the site to the reset list
        pf_resetlist[pf_nresetlist]=i;
        pf_resetlistvals[pf_nresetlist]=val;
        
        pf_nresetlist++;
      }
    iarg += 2;
  }
    
  //make sure no site is reset more than once to a different value
  for (int i=0;i<pf_nresetlist;i++)
    for (int j=0;j<pf_nresetlist;j++) {
      if (i==j) continue;
      else if (pf_resetlist[i] == pf_resetlist[j] && 
               pf_resetlistvals[i] != pf_resetlistvals[j])
        error->all(FLERR,"Invalid reset_phasefield command: attempting to "
                   "reset a site more than once to different values");
    }
  
  //arguments have been parsed, set flag
  pf_resetfield = true;
  
  // print statistics
  
  tagint nbig = pf_nresetlist;
  tagint allcount;
  MPI_Allreduce(&nbig,&allcount,1,MPI_SPK_TAGINT,MPI_SUM,world);
  
  if (domain->me == 0) {
    if (screen) fprintf(screen,"Setting reset phasefield list ...\n"
                        "  " TAGINT_FORMAT " sites identified to reset\n",
                        allcount);
    if (logfile) fprintf(logfile,"Setting reset phasefield list ...\n"
                         "  " TAGINT_FORMAT " sites identified to reset\n",
                         allcount);
  }
}

/* ---------------------------------------------------------------------- */
//binary search to find nearest value. code adapted from
// http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html
int AppPottsPhaseFieldEric::binarySearch(double *sortedArray, 
                                     int first, int last, double key) {
  //this returns the index for the nearest value
  while (1) {
    int mid = (last - first) / 2;  // compute mid point.
    if (mid <= 0) {
      if (abs(key - sortedArray[first]) <= abs(key - sortedArray[last]))
        return first;
      else
        return last;
    }
    mid = first + mid;
    if (key > sortedArray[mid]) 
      first = mid;  // repeat search in top half.
    else if (key < sortedArray[mid]) 
      last = mid; // repeat search in bottom half.
    else
      return mid;     // found it. return position
  }
}

/* ----------------------------------------------------------------------
 get app properties or in this case the pointer to the table
 ------------------------------------------------------------------------- */

void *AppPottsPhaseFieldEric::extract_app(char * name) {
  
  int i=find_phase(name);
  if (i < 0)
    error->all(FLERR,"Cannot find the phase");
  else 
    return (void *) phases[i]->table;
  
  return (void *) NULL;
}
