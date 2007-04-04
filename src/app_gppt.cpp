/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_gppt.h"
#include "solve_gppt.h"
#include "spk.h"
#include "solve.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include <iostream>
#include "node.h"
#include "plus_node.h"
#include "minus_node.h"
#include "times_node.h"
#include "divide_node.h"
#include "const_node.h"
#include "var_node.h"
#include "pow_node.h"
#include "tree.h"
#include "population.h"

using namespace std;
using namespace SPPARKS;


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */
AppGppt::AppGppt(SPK *spk, int narg, char **arg) : App(spk, narg, arg)
{
  if (narg != 2) error->all("Invalid app_style GPPT command");

  if (screen) fprintf(screen,"GPPT simulation.\n");
  if (logfile)fprintf(logfile,"GPPT simulation.\n");

  // define proc layout
  
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  procup = me + 1;
  if (procup > nprocs - 1) procup = 0;
  procdown = me - 1;
  if (procdown < 0) procdown = nprocs - 1;

  ntimestep = 0;
  time = 0.0;
  stoptime = 0.0;
  // classes needed by this app
  seed = atoi(arg[1]);

  if (random != NULL)delete random;
  random = new RandomPark(seed);

  if (timer != NULL)delete timer;
  timer = new Timer(spk);

  if (tree != NULL)delete tree;
  tree = new Tree();

  if (pop != NULL)delete pop;
  //pop = new *Population();

  if (fit != NULL)delete fit;
  fit = new Fitness(spk);

  stats_delta = 0;
  max_tree_depth = 1;
  pop_size = 0;
  pop_num = 0;
  num_variables = 0;

  best_fitness = 1.0e12;
  best_pop = 0;

  thi = 0.0;
  tlo = 0.0;
  swap_freq = 0;

  temp_init = "linear";
}
/* ---------------------------------------------------------------------- */
AppGppt::~AppGppt()
{
  //  best_tree->clear();
  delete random;
  delete tree;
  delete fit;
  delete [] fitness;
  delete [] pop;
  delete [] var;
  for (int b = 0; b < buf_size; b++) delete [] buffer[b];
  delete [] buffer;
  delete [] tight_buf;
}
/* ---------------------------------------------------------------------- */
void AppGppt::init()
{
  Node *current_root;
  double current_t;
  int pop_per_proc;

  if(num_variables > 0 & range > 0)
    tree->init(seed, const_lo, const_hi, num_variables, max_tree_depth);
  else {
    error->all("Number of variables or constant range not set.");
    return;
  }

  buf_size = 2;
  word_size = 16;
  for (int i = 0; i < max_tree_depth; i++) buf_size *= 2;
  buffer = new char*[buf_size];
  for (int b = 0; b < buf_size; b++) buffer[b] = new char[word_size];
  tight_buf = new char[buf_size*word_size];

  if(num_variables == fit->get_num_var()) fit->init();
  else
    error->all("Number of variables inconsistent in tree and fitness.");

  if(swap_freq == 0) error->all("Swap frequency is not set.");



  //distribute populations to procs

  if(pop_size > 0 & pop_num > 0){

    pop_per_proc = pop_num/nprocs;
    if(pop_num > pop_per_proc * nprocs) pop_per_proc++;
    if(pop_per_proc * nprocs > pop_num)
      if(me == nprocs - 1) pop_per_proc = pop_num - pop_per_proc*nprocs;
    pop_num = pop_per_proc;

    pop = new Population[pop_num];
    current_t = thi;
    for(int p = 0; p < pop_num; p++){
      pop[pop_num-p-1].init(seed, pop_size, 
			  tree, fit, max_tree_depth, current_t);
      if (strcmp(temp_init.c_str(),"log")==0) current_t = (tlo + current_t)/2; 
      else if (strcmp(temp_init.c_str(),"linear")==0)
	current_t -= (thi-tlo)/static_cast<double>(pop_num);
    }
  }
  else {
    error->all("Population parameters not set.");
    return;
  }

  get_best_fitness();
  
  cout <<"TREE TEST: "<<endl;
//   double test_var[1];
//   test_var[0] = 2.5;
//   best_tree->write(screen);
//   cout <<endl;
//   best_tree->write_tex(screen);
//   cout <<endl;
//   write_tex_header(best_tree);
//   cout <<endl;
//   cout <<"Fitness: "<<best_tree->go(test_var)<<endl;
  
  MinusNode *root = new MinusNode();
  PowNode *twelve = new PowNode();
  PowNode *six = new PowNode();
  //   DivideNode *divleft = new DivideNode();
  ConstNode *exp_left = new ConstNode();
  //   DivideNode *divright = new DivideNode();
  ConstNode *exp_right = new ConstNode();
  VarNode *var_left = new VarNode();
  VarNode *var_right = new VarNode();
  
  root->set_left_child(static_cast<Node *>(twelve));
  root->set_right_child(static_cast<Node *>(six));
  
  twelve->set_left_child(static_cast<Node *>(var_left));
  twelve->set_right_child(static_cast<Node *>(exp_left));
  
  var_left->set_value(0);
  exp_left->set_value(-12.0);
  
  six->set_left_child(static_cast<Node *>(var_right));
  six->set_right_child(static_cast<Node *>(exp_right));
  
  var_right->set_value(0);
  exp_right->set_value(-6.0);
  
  
  write_tex_header(root);
  cout <<endl;
  cout << "Fitness of the test tree = "<< fit->compute(root)<<endl;
  
  root->clear();
  
//   cout <<"R**R TEST: "<<endl;
//   PowNode *proot = new PowNode();
//   VarNode *pvar_left = new VarNode();
//   VarNode *pvar_right = new VarNode();
//   proot->set_left_child(static_cast<Node *>(pvar_left));
//   proot->set_right_child(static_cast<Node *>(pvar_right));
//   pvar_left->set_value(0);
//   pvar_right->set_value(0);

  cout <<"CONST TEST: "<<endl;
  
  ConstNode *proot = new ConstNode();
  proot->set_value(0.98274);
  
  write_tex_header(proot);
  cout <<endl;
  cout << "Fitness of the test tree = "<< fit->compute(proot)<<endl;
  
  proot->clear();
  
  if (screen) {
    fprintf
      (screen,"Step Time BestPop BestFitness const var oper branch cross swap");
    fprintf(screen,"\n");
  }
  if (logfile) {
    fprintf
      (logfile,"Step Time BestPop BestFitness const var oper branch cross swap");
    fprintf(logfile,"\n");
  }
  stats();
  
  // setup future calls to stats()
  
  stats_time = time + stats_delta;
  if (stats_delta == 0.0) stats_time = stoptime;
}
/* ---------------------------------------------------------------------- */
void AppGppt::write_tex_header(Node *tree)
{
  if (screen){
    fprintf(screen,"######################## Tree ############.\n");
    fprintf(screen,"\\documentclass[letter,10pt]{article} \n");
    fprintf(screen,"\\begin{document}\n");
    fprintf(screen,"$$\n");
    tree->write_tex(screen);
    fprintf(screen,"\n");
    fprintf(screen,"$$\n");
    fprintf(screen,"\\end{document}\n");
  }
}
/* ---------------------------------------------------------------------- */
void AppGppt::write_tree(int t)
{
  fprintf(screen,"######################## Tree %d############.\n",t);
  trees[t]->write(screen);
  fprintf(screen,"\n");
}

/* ---------------------------------------------------------------------- */
void AppGppt::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  else if (strcmp(command,"tree_type") == 0) set_tree_type(narg,arg);
  else if (strcmp(command,"population") == 0) set_population(narg,arg);
  else if (strcmp(command,"tempering") == 0) set_tempering(narg,arg);
  else if (strcmp(command,"fitness") == 0) set_fitness(narg,arg);
  else if (strcmp(command,"run") == 0) run(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);

  else error->all("Invalid command");
}
/* ----------------------------------------------------------------------
   perform a run
------------------------------------------------------------------------- */
void AppGppt::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

  // error check

  if (solve == NULL) error->all("No solver class defined");
  //  if (tree == NULL) error->all("No tree class defined");
  //  if (fit == NULL) error->all("No fitness class defined");

  // init classes used by this app
  
  init();
  //  solve->init(nevents,propensity);
  timer->init();

  // perform the run

  iterate();

  // final statistics

  Finish finish(spk);
}
/* ----------------------------------------------------------------------
   iterate on solver
------------------------------------------------------------------------- */
void AppGppt::iterate()
{
  int d,ievent;
  double dt = 1.0;
  int p;
  int done = 0;
  int par = 0;

  timer->barrier_start(TIME_LOOP);

  while (!done) {

    ntimestep++;
    timer->stamp();
    //    ievent = solve->event(&dt);
    for(p = 0; p < pop_num; p++) pop[p].build_new_population();

    timer->stamp(TIME_SOLVE);

    // update propensity table
    // inform solver of changes

    //update event propensity

    //    solve->update(ievent, propensity);

    //update mutation weights
    for(p = 0; p < pop_num; p++) pop[p].update_weights();
    //swap trees
    if (static_cast<int>(time)%swap_freq == 0) {
      swap_trees(par);
      par += 1;
      par = par%2;
    }

    if (static_cast<int>(time)%ndump == 0) dump(); 
    //update dependencies
    
    //    solve->update(ndepends[ievent],depends[ievent],propensity);
    
    timer->stamp(TIME_UPDATE);
    // update time by dt

    time += dt;
    if (time >= stoptime) done = 1;
    else if (ievent < 0) done = 1;

    if (time > stats_time || done) {
      timer->stamp();
      stats();
      stats_time += stats_delta;
      timer->stamp(TIME_OUTPUT);
    }
  }

  timer->barrier_stop(TIME_LOOP);
}
/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */
void AppGppt::stats()
{
  int i;
  double *weight;
  int nw = 0;

  weight = pop[best_pop].get_weights(nw);

  if (screen) {
    fprintf(screen,"%d %g %d %g ",ntimestep,time, best_pop, 
	    get_best_fitness());
    fprintf(screen,"%g %g %g %g %g %g ",weight[0],weight[1], weight[2],
	    weight[3], weight[4], weight[5]);
    fprintf(screen,"\n \n");
  }
  if (logfile) {
    fprintf(logfile,"%d %g %d %g  ",ntimestep,time, best_pop,
	    get_best_fitness());
    fprintf(logfile,"%g %g %g %g %g %g ",weight[0],weight[1], weight[2],
	    weight[3], weight[4], weight[5]);
    fprintf(logfile,"\n \n");
  }

  if (screen){
    best_tree->write(screen);
    fprintf(screen, "\n \n");
    best_tree->write_tex(screen);
    fprintf(screen, "\n \n");
  }
  if (logfile){
    best_tree->write(logfile);
    fprintf(logfile, "\n \n");
    best_tree->write_tex(logfile);
    fprintf(logfile, "\n \n");
  }
}
/* ----------------------------------------------------------------------
   proc 0 writes dump header
------------------------------------------------------------------------- */

void AppGppt::dump_header()
{
  if(me == 0){
    if (fp) {
      fprintf
	(fp,"Step Time Temp1 Mean1 Stdev1 Best1 Temp2 Mean2 Stdev2 Best2 ");
      fprintf(fp,"\n");
    }

  }
}
/* ----------------------------------------------------------------------
   dump a snapshot of two populations stats
------------------------------------------------------------------------- */

void AppGppt::dump()
{
  double mean1, mean2;
  double stdev1, stdev2;
  double best1, best2;
  double temp1, temp2;

  mean1 = pop[dpop1].get_stats(stdev1);
  mean2 = pop[dpop2].get_stats(stdev2);
  best1 = pop[dpop1].get_best_fitness();
  best2 = pop[dpop2].get_best_fitness();
  temp1 = pop[dpop1].get_temperature();
  temp2 = pop[dpop2].get_temperature();

  if (me == 0) {
    fprintf(fp,"%d %g ",ntimestep,time);
    fprintf(fp,"%g %g %g %g %g %g %g %g ",temp1,mean1,stdev1,best1,
	    temp2,mean2,stdev2,best2);
    fprintf(fp,"\n");
  }

}
/* ---------------------------------------------------------------------- */
void AppGppt::set_tree_type(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal tree_type command");
  max_tree_depth = atoi(arg[0]);

  const_lo = atof(arg[1]);
  const_hi = atof(arg[2]);
  range = const_hi - const_lo + 1;
  num_variables = atoi(arg[3]);
  if(screen){
    fprintf(screen,"Tree depth = %d\n", max_tree_depth);
    fprintf(screen,"Constants bounded between %g and %g \n",const_lo, const_hi);
    fprintf(screen,"Sensor variables = %d\n", num_variables);
  }
  if(logfile){
    fprintf(logfile,"Tree depth = %d\n", max_tree_depth);
    fprintf(logfile,"Constants bounded between %g and %g \n",const_lo, const_hi);
    fprintf(logfile,"Sensor variables = %d\n", num_variables);
  }
}
/* ---------------------------------------------------------------------- */
void AppGppt::set_population(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal population command");
  pop_num = atoi(arg[0]);
  pop_size = atoi(arg[1]);
  if(screen){
    fprintf(screen,"Population number %d\n", pop_num);
    fprintf(screen,"Population size %d\n", pop_size);
  }
  if(logfile){
    fprintf(logfile,"Population number %d\n", pop_num);
    fprintf(logfile,"Population size %d\n", pop_size);
  }
}
/* ---------------------------------------------------------------------- */
void AppGppt::set_tempering(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal tempering command");

  temp_init.assign(arg[0]);

  if(!(temp_init == "log") & !(temp_init == "linear"))
    error->all("Illegal temperature distribution");

  tlo = atof(arg[1]);
  thi = atof(arg[2]);
  swap_freq = atoi(arg[3]);

  if(screen){
    fprintf(screen,"Setting initial T distribution type: ");
    fprintf(screen,temp_init.c_str());
    fprintf(screen,"\n");
    fprintf(screen,"Lower bound T = %g\n", tlo);
    fprintf(screen,"Upper bound T = %g\n", thi);
    fprintf(screen,"Swap frequency = %d\n", swap_freq);
  }
  if(logfile){
    fprintf(logfile,"Setting initial T distribution type:\n");
    fprintf(logfile,temp_init.c_str());
    fprintf(logfile,"\n");
    fprintf(logfile,"Lower bound T = %g\n", tlo);
    fprintf(logfile,"Upper bound T = %g\n", thi);
    fprintf(logfile,"Swap frequency = %d\n", swap_freq);
  }
}
/* ---------------------------------------------------------------------- */
void AppGppt::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);

  if(screen){
    fprintf(screen,"Stats frequency = %g\n", stats_delta);
  }
  if(logfile){
    fprintf(logfile,"Stats frequency = %g\n", stats_delta);
  }
}
/* ---------------------------------------------------------------------- */

void AppGppt::set_dump(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal dump command");
  ndump = atoi(arg[0]);
  if (me == 0) {
    fp = fopen(arg[1],"w");
    if (fp == NULL) error->one("Cannot open dump file");
  }
  //two dump populations index
  dpop1 = atoi(arg[2]);
  dpop2 = atoi(arg[3]);
}
/* ---------------------------------------------------------------------- */
void AppGppt::set_fitness(int narg, char **arg)
{
  fit->set_fitness(narg, arg, me);
}
/* ---------------------------------------------------------------------- */
double AppGppt::get_best_fitness()
{
  double current_best_fitness;

  for(int p = 0; p < pop_num; p++) {
    current_best_fitness = pop[p].get_best_fitness();
    if(current_best_fitness < best_fitness){
      best_fitness = current_best_fitness;
      best_tree = pop[p].get_best_tree();
      best_pop = p;
    }
  }
  return best_fitness;
}

/* ---------------------------------------------------------------------- */
void AppGppt::swap_trees(int par)
{
  Population *curr_pop;
  Population *curr_nbr;
  int in_size = 0;
  char *inc_buf;
  int s;
  int nbr;
  double nbr_temp;

  //populations pack buffers with best trees
  //send to mailbox
  //neigbors receive
  for(int p = 0; p < pop_num; p++) {
    //pack buffers for send
    curr_pop = &pop[p];
    curr_pop->pack_buffer(curr_pop->get_best_tree());
    //receive at mailbox
    inc_buf = curr_pop->get_buffer(in_size);
    s=0;
    while(1==1){
      tight_buf[s] = inc_buf[s];
      if(inc_buf[s] == '#') break;
      s++;
    }
    //determine neighbor population index
    if (par == 0) nbr = (p+1)%pop_num;
    else if (par == 1) {
      if(p > 0) nbr = p-1;
      else nbr = pop_num - 1; 
    }
    curr_nbr = &pop[nbr];
    //unpack buffers upon receipt
    curr_nbr->unpack_buffer(tight_buf);
    //create tree from buffer and check acceptance
    nbr_temp = curr_pop->get_temperature();
    curr_nbr->swap_nbr(par, nbr_temp);
  }
}
/* ---------------------------------------------------------------------- */
