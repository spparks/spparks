/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_custom.h"
#include "flip_event.h"
#include "swap_event.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppCustom::AppCustom(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // parse arguments

  if (narg < 3) error->all("Invalid app_style custom command");

  seed = atoi(arg[1]);
  //  nspins = atoi(arg[2]);


  random = new RandomPark(seed);

  options(narg-2,&arg[2]);

  // define lattice and partition it across processors
  
  create_lattice();
  sites = new int[1 + 10*maxneigh];

  // assign variable names

  if (ninteger == 0 & ndouble == 0)
    error->all("Invalid site specification in app_style custom");

  temperature = -1;
  energyxp = NULL;

  //create a "zero" tree for initialization

  char **zargs;
  zargs = new char*[2];
  zargs[0] = new char[5];
  strcpy(zargs[0],"ZERO");
  zargs[1] = new char[5];
  strcpy(zargs[1],"0.0");
  set_expression(2, zargs);

  //init arrays and vars for the state name lookup

  n_int_var = 0;
  n_dbl_var = 0;
  name_int = new char*[ninteger];
  name_dbl = new char*[ndouble];

  //init arrays for nbrhood union for updates

  update_numneigh = new int[nlocal];
  update_neighbor = new int*[nlocal];
  for (int i = 0; i < nlocal; i++){
    update_numneigh[i] = 0;
    update_neighbor[i] = NULL;
  }

  tr = NULL;

  // initialize hash for n-proc independent random assignment

  for (int i = 0; i < nlocal; i++)
    hash.insert(std::pair<int,int> (id[i],i));
}

/* ---------------------------------------------------------------------- */

AppCustom::~AppCustom()
{
  delete random;
  delete [] sites;
  delete comm;
  delete [] name_int;
  delete [] name_dbl;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppCustom::site_energy(int i)
{
  return energyxp->sum(i);
}

/* ----------------------------------------------------------------------
   randomly pick new state for site
------------------------------------------------------------------------- */

void AppCustom::site_pick_random(int i, double ran)
{
//   int iran = (int) (nspins*ran) + 1;
//   if (iran > nspins) iran = nspins;
//   return iran;
}

/* ----------------------------------------------------------------------
   randomly pick new state for site from neighbor values
------------------------------------------------------------------------- */

void AppCustom::site_pick_local(int i, double ran)
{
//   int iran = (int) (numneigh[i]*ran) + 1;
//   if (iran > numneigh[i]) iran = numneigh[i];
//   return spin[neighbor[i][iran]];
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site
   based on einitial,efinal for each possible event
   if no energy change, propensity = 1
   if downhill energy change, propensity = 1
   if uphill energy change, propensity set via Boltzmann factor
   if proc owns full domain, there are no ghosts, so ignore full flag
------------------------------------------------------------------------- */

double AppCustom::site_propensity(int i, int full)
{
  double prop_sum = 0;
  int nevent = events.size();

  for (int ev = 0; ev < nevent; ev++)
    prop_sum += events[ev]->site_propensity(i);

  return prop_sum;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
   update propensities of all affected sites
   if proc owns full domain, there are no ghosts, so ignore full flag
   if proc owns sector, ignore neighbor sites that are ghosts
------------------------------------------------------------------------- */

void AppCustom::site_event(int i, int full)
{
  int j,k,m,isite,value;
  double efinal;
  int ev;

  // pick one event from total propensity

  double threshold = random->uniform() * propensity[i2site[i]];

  int nevent = events.size();

  int second_site = -1;
  for (ev = 0; ev < nevent; ev++){
    int val = events[ev]->site_event(i, threshold);
    if(val == 1) {
      break;
    }
    if(val > 99){
      second_site = val - 100;
      break;
    }
  }

  // compute propensity changes for self and neighbor sites

  int nsites = 0;
  isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i,full);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m,full);
  }

  for (j = 0; j < update_numneigh[i]; j++) {
    m = update_neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m,full);
  }

  //additional nbrs for a two-site move (like swap)

  if(second_site > -1){
    bool included = false;
    int ssite;
    int s;
    ssite = i2site[second_site];
    for (s=0;s<nsites;s++) if(ssite == sites[s]) included = true;
    if(!included) {
      sites[nsites++] = ssite;
      propensity[ssite] = site_propensity(i,full);
      included = false;
    }
    for (j = 0; j < numneigh[i]; j++) {
      m = neighbor[second_site][j];
      ssite = i2site[m];
      if (ssite < 0) continue;
      for (s=0;s<nsites;s++) if(ssite == sites[s]) included = true;
      if(!included) {
	sites[nsites++] = ssite;
	propensity[ssite] = site_propensity(m,full);
	included = false;
      }
    }
    
    for (j = 0; j < update_numneigh[i]; j++) {
      m = update_neighbor[second_site][j];
      ssite = i2site[m];
      if (ssite < 0) continue;
      for (s=0;s<nsites;s++) if(ssite == sites[s]) included = true;
      if(!included) {
	sites[nsites++] = ssite;
	propensity[ssite] = site_propensity(m,full);
	included = false;
      } 
    }
  }
  
  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
  clear mask values of site and its neighbors
  OK to clear ghost site mask values
------------------------------------------------------------------------- */

void AppCustom::site_clear_mask(char *mask, int i)
{
  mask[i] = 0;
  for (int j = 0; j < numneigh[i]; j++) mask[neighbor[i][j]] = 0;
}

/* ----------------------------------------------------------------------
   parse commands
------------------------------------------------------------------------- */

void AppCustom::input_app(char *command, int narg, char **arg)
{
  //  fprintf(screen,"command = %s \n",command);
  if      (strcmp(command,"var")  == 0) set_variable(narg,arg);
  else if (strcmp(command,"init") == 0) init_variable(narg,arg);
  else if (strcmp(command,"expr") == 0) set_expression(narg,arg);
  else if (strcmp(command,"neighborhood") == 0) set_neighborhood(narg,arg);
  else if (strcmp(command,"event") == 0) set_event(narg,arg);
  else if (strcmp(command,"propensity") == 0) set_propensity(narg,arg);
  else if (strcmp(command,"energy") == 0) set_energy(narg,arg);
  else error->all("Command not recognized by this application");
}

/* ----------------------------------------------------------------------
   add variable
------------------------------------------------------------------------- */

void AppCustom::set_variable(int narg, char **arg)
{
  if (narg != 2) error->one("Illegal variable command");
  int type = 0;

  if      (strcmp(arg[1],"integer") == 0) type = 1;
  else if (strcmp(arg[1],"double") == 0 ) type = 2;
  else error->one("Unknown variable type.");

  if(type == 1){
    if(n_int_var+1 > ninteger) 
      error->all("Not enough integer storage on site.");
    name_int[n_int_var] = new char[strlen(arg[0])+1];
    strcpy(name_int[n_int_var], arg[0]);
    n_int_var ++;
  }
  else if(type == 2){
    if(n_dbl_var+1 > ndouble) 
      error->all("Not enough double storage on site.");
    name_dbl[n_dbl_var] = new char[strlen(arg[0])+1];
    strcpy(name_dbl[n_dbl_var], arg[0]);
    n_dbl_var ++;
  }

}

/* ----------------------------------------------------------------------
   initialize variable
------------------------------------------------------------------------- */

void AppCustom::init_variable(int narg, char **arg)
{
  if (narg < 2) error->one("Illegal init command");
  int type = 0;

  int np = get_variable(arg[0],type);
  //  fprintf(screen, "Init variable %s of type %d to %s \n",arg[0],type,arg[2]);

  double lo, hi;
  double cd = 0.0;
  int ci = 0;
  int v;

  if(strcmp(arg[1],"constant")==0){
    if(type == 2){
      cd = atof(arg[2]);
      for(v = 0; v < nlocal; v++) darray[np][v] = cd;
    }
    else if(type == 1) {
      ci = atoi(arg[2]);
      for(v = 0; v < nlocal; v++) iarray[np][v] = ci;
    }
  }
  else if(strcmp(arg[1],"linear")==0){
    lo =  atof(arg[2]);
    hi =  atof(arg[3]);
  }
  else if(strcmp(arg[1],"random")==0){

    lo =  atof(arg[2]);
    hi =  atof(arg[3]);

    int isite;
    double dsite;

    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      if     (type == 1) {
	int range = (int)hi - (int)lo + 1;
	isite = (int)lo - 1 + random->irandom(range);
      }
      else if(type == 2) {
	dsite = lo + (hi-lo)*random->uniform();
      }
      loc = hash.find(iglobal);
      if (loc != hash.end()) {
	if     (type == 2) darray[np][loc->second] = dsite;
	else if(type == 1) iarray[np][loc->second] = isite;
      }
    }
  }
  else error->all("Unknown distribution type in init.");
//   if(type == 1)
//     for(v = 0; v < nlocal; v++) 
//       fprintf(screen,"int %d = %d\n", v, iarray[np][v]);
//   if(type == 2)
//     for(v = 0; v < nlocal; v++) 
//       fprintf(screen,"dbl %d = %g\n", v, darray[np][v]);
}

/* ----------------------------------------------------------------------
   find variable by name
------------------------------------------------------------------------- */

int AppCustom::get_variable(char *name, int &type)
{
  for (int i = 0; i < n_int_var; i++){
    if(strcmp(name, name_int[i])==0) {
      type = 1; 
      return i;
    }
  }
  for (int i = 0; i < n_dbl_var; i++){
    if(strcmp(name, name_dbl[i])==0) {
      type = 2; 
      return i;
    }
  }
  error->all("Variable not found.");
  return -1;
}

/* ----------------------------------------------------------------------
   add expression
------------------------------------------------------------------------- */

void AppCustom::set_expression(int narg, char **arg)
{
  if (narg < 2) error->one("Illegal expression command");

  if (tr == NULL) {
    tr = new Tree();
    tr->init(0.1, 1, 1, 10);
    tr->set_lattice(iarray, darray);
    tr->set_names(name_int, n_int_var, name_dbl, n_dbl_var);
  }
  Expression *xpr_local = new Expression();
  xpr_local->set_tree(tr);
  xpr_local->set_expression(narg, arg);
  expr.push_back(xpr_local);

}

/* ----------------------------------------------------------------------
   add neighborhood
------------------------------------------------------------------------- */

void AppCustom::set_neighborhood(int narg, char **arg)
{
  int test_nb=-1;
  int nb = -1;

  if (narg < 2) error->one("Illegal neighborhood command");

  Neighborhood *nbhd_local = 
    new Neighborhood(arg[0], numneigh, neighbor, nlocal, maxneigh);

  nbhd_local->init(narg, arg);
  nbhd.push_back(nbhd_local);

  //add to the union of update neighbors
  for(int st=0;st<nlocal;st++){
    for(int nbn = 0;nbn < nbhd_local->numneigh[st]; nbn ++){
      test_nb = nbhd_local->neighbor[st][nbn];
      for(nb = 0; nb < numneigh[st]; nb++){
	if(test_nb == neighbor[st][nb])break;
      }
      if(nb == numneigh[st]){
	update_numneigh[st]++;
	update_neighbor[st] = (int *)
	  realloc((char *)update_neighbor[st],
		  update_numneigh[st] * sizeof(int));
	update_neighbor[st][update_numneigh[st]-1] = test_nb;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   find a neighborhood
------------------------------------------------------------------------- */

Neighborhood *AppCustom::get_neighborhood(char *name_in)
{
  for (int nb = 0; nb < nbhd.size(); nb++){
    if(strcmp(name_in, nbhd[nb]->get_name())==0) return nbhd[nb];
  }
  error->all("Neighborhood not found.");
  return NULL;
}

/* ----------------------------------------------------------------------
   find an expression
------------------------------------------------------------------------- */

Expression *AppCustom::get_expression(char *name_in)
{
  for (int xp = 0; xp < expr.size(); xp++){
    if(strcmp(name_in, expr[xp]->get_name())==0) return expr[xp];
  }
  error->all("Expression not found.");
  return NULL;
}

/* ----------------------------------------------------------------------
   add an event
------------------------------------------------------------------------- */

void AppCustom::set_event(int narg, char **arg)
{
  Event *current_event;
  int vtype = 0;
  int numvar = 0;
  int nword = 0;

  if (narg < 3) error->one("Illegal event command");

  if(temperature == -1) 
    error->one("Need to set temeprature before the event");

  //parse event description here
  //event_name event_type
  if(strcmp(arg[1],"flip")==0) {
    current_event = new FlipEvent();
    events.push_back(current_event);
  }
  else if(strcmp(arg[1],"swap")==0) {
    current_event = new SwapEvent();
    events.push_back(current_event);
  }
  current_event->set_name(arg[0]);
  //variable name
  numvar = get_variable(arg[2],vtype);
  if (vtype==1)
    current_event->set_variable(vtype, iarray[numvar]);
  else if (vtype==2)
    current_event->set_variable(vtype, darray[numvar]);
  else if(vtype==-1)
    error->all("Event variable not found.");

  //consistency checking
  if(strcmp(arg[1],"flip")==0 & vtype !=1)
    error->all("Flip event operates on integer variables.");
  if(strcmp(arg[1],"swap")==0 & vtype !=1)
    error->all("Swap event operates on integer variables.");

  nword = 3; //keep track of optional arguments

  //update source
  //to_nbr neighborhood_name
  //enum hi lo
  if(strcmp(arg[nword],"to_nbr")==0){
    nword++;
    current_event->set_neighborhood(get_neighborhood(arg[nword]));
    current_event->set_src(10);
    nword++;
  }
  else if (strcmp(arg[nword],"enum")==0){
    if(strcmp(arg[1],"swap")==0)
      error->all("Swap event cannot have enum source.");
    nword++;
    current_event->set_enum_bounds(atoi(arg[nword]),atoi(arg[nword+1]));
    current_event->set_src(11);
    nword++; nword++;
  }
  else if(strcmp(arg[1],"swap")==0)
      error->all("Swap event must have to_nbr source.");

  current_event->set_temperature(temperature);
  //default propensity is a ZERO tree
  Expression *my_expression = get_expression("ZERO");
  current_event->set_propensity(1,my_expression);

  current_event->init();
}

/* ----------------------------------------------------------------------
   find an event
------------------------------------------------------------------------- */

Event *AppCustom::get_event(char *name_in)
{
  for (int ev = 0; ev < events.size(); ev++){
    if(strcmp(name_in, events[ev]->get_name())==0) return events[ev];
  }
  error->all("Event not found.");
  return NULL;
}

/* ----------------------------------------------------------------------
   set propensity for an event
------------------------------------------------------------------------- */

void AppCustom::set_propensity(int narg, char **arg)
{
  int ptype = -1;

  if (narg < 2) error->one("Illegal propensity command");
  Event *my_event = get_event(arg[0]);
  Expression *my_expression = get_expression(arg[2]);

  if(strcmp(arg[1],"boltz")==0) ptype = PROP_BOLTZ;
  else if(strcmp(arg[1],"eval")==0) ptype = PROP_EVAL;
  else if(strcmp(arg[1],"nb_sum")==0) ptype = PROP_NB_SUM;

  my_event->set_propensity(ptype, my_expression);
}

/* ----------------------------------------------------------------------
   set energy expression
------------------------------------------------------------------------- */

void AppCustom::set_energy(int narg, char **arg)
{
  if (narg < 1) error->one("Illegal energy command");
  energyxp = get_expression(arg[0]);
}

