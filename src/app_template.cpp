/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_template.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include <iostream>
#include "tree.h"
#include "output.h"

using namespace SPPARKS;
using namespace std;



// #define MIN(a,b) ((a) < (b) ? (a) : (b))
// #define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppTemplate::AppTemplate(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
{
  if (narg != 1) error->all("Invalid app_style template command");

  // default settings

  ntimestep = 0;
  time = 0.0;
 
  st = NULL;
  tr = NULL;
  xpr = NULL;
  nxpr = 0;

}

/* ---------------------------------------------------------------------- */

AppTemplate::~AppTemplate()
{
}

/* ---------------------------------------------------------------------- */

void AppTemplate::init()
{


  //check if set state node attributes

  //check if set edge sets 

  //check if set expressions

  //check if set event descriptors

  //use event descriptors to generate initial event table

  // Initialize output

  output->init(time);

}

/* ---------------------------------------------------------------------- */

void AppTemplate::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");

  if (strcmp(command,"run") == 0) run(narg,arg);

  else if (strcmp(command,"stats") == 0) output->set_stats(narg,arg);

  else if (strcmp(command,"state") == 0) set_state(narg,arg);

  else if (strcmp(command,"var") == 0) set_variable(narg,arg);

  else if (strcmp(command,"init") == 0) init_variable(narg,arg);

  else if (strcmp(command,"expr") == 0) set_expression(narg,arg);


  else error->all("Invalid command");
}
/* ----------------------------------------------------------------------
   create state
------------------------------------------------------------------------- */
void AppTemplate::set_state(int narg, char **arg)
{
  if (narg != 2) error->one("Illegal state command");
  st = new State(atoi(arg[0]),atoi(arg[1]));
  tr = new Tree();
  tr->init(0.1, 1, 1, 10);
  tr->set_state(st);

}
/* ----------------------------------------------------------------------
   add variable
------------------------------------------------------------------------- */
void AppTemplate::set_variable(int narg, char **arg)
{
  if (narg != 2) error->one("Illegal variable command");
  if (st == NULL) error->one("State must be defined before variables.");

  int type = 0;

  if      (strcmp(arg[1],"integer") == 0) type = 1;
  else if (strcmp(arg[1],"double") == 0 ) type = 2;
  else error->one("Unknown variable type.");
  st->add_attribute(arg[0], type);
}
/* ----------------------------------------------------------------------
   add variable
------------------------------------------------------------------------- */
void AppTemplate::init_variable(int narg, char **arg)
{
  if (narg < 2) error->one("Illegal init command");
  st->init_var(narg, arg);
  
}
/* ----------------------------------------------------------------------
   add expression
------------------------------------------------------------------------- */
void AppTemplate::set_expression(int narg, char **arg)
{
  if (narg < 2) error->one("Illegal init command");
  if (st == NULL) error->one("State must be defined before expressions.");

  Xpression *xpr_local = new Xpression();
  xpr_local->set_state(st);
  xpr_local->set_tree(tr);
  xpr_local->set_expression(narg, arg);
  expr.push_back(xpr_local);
  nxpr++;

}
/* ----------------------------------------------------------------------
   perform a run
------------------------------------------------------------------------- */

void AppTemplate::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

  for (int e = 0; e < nxpr; e++ ){
    cout <<"next expression: "<<endl;
    for(int site = 0; site < st->get_size(); site++)
      cout << expr[e]->eval(site)<<endl;
  }

  //%%%%%%%%%% diagnostic for state content
  //   char *tmp;
  //   strcpy(tmp,"s2");
  //   int type;
  
  //   int *tempi = (int *)st->get_attribute(tmp, type);
  
  //   cout <<"testing contents:"<<endl;
  //   for(int i = 0; i < 100; i++) cout << i<<"  "<<tempi[i]<<endl;
  //%%%%%%%%%% end diagnostic

  // error check

  if (solve == NULL) error->all("No solver class defined");

  // init classes used by this app
  
  init();
  //  solve->init(nreactions,propensity);
  timer->init();

  // perform the run

  iterate();

  // final statistics

  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on Gillespie solver
------------------------------------------------------------------------- */

void AppTemplate::iterate()
{
  int m,ireaction;
  double dt;

  int done = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {
    ntimestep++;

    timer->stamp();
    //    ireaction = solve->event(&dt);
    timer->stamp(TIME_SOLVE);

    //    solve->update(ndepends[ireaction],depends[ireaction],propensity);
 
    time += dt;
    if (time >= stoptime) done = 1;

    // Do output

    output->compute(time,done);
    timer->stamp(TIME_OUTPUT);

  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppTemplate::stats(char *strtmp)
{
  sprintf(strtmp," %10d %10g",ntimestep,time);
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppTemplate::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s","Step","Time");
}

/* ---------------------------------------------------------------------- */

void AppTemplate::set_stats(int narg, char **arg)
{
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"your_option_here") == 0) {
      iarg++;
      if (iarg < narg) {
	int itmp = atoi(arg[iarg]);
      } else {
	error->all("Illegal stats command");
      }
    }
    iarg++;
  }
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
