/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "app_gppt.h"
#include "solve_gppt.h"
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

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppGPPT::AppGPPT(SPK *spk, int narg, char **arg) : App(spk, narg, arg)
{
  if (narg != 2) error->all("Invalid app_style GPPT command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (me == 0) {
    if (screen) fprintf(screen,"GPPT simulation.\n");
    if (logfile)fprintf(logfile,"GPPT simulation.\n");
  }

  // define proc layout
  
  procup = me + 1;
  if (procup > nprocs - 1) procup = MPI_PROC_NULL;
  procdown = me - 1;
  if (procdown < 0) procdown = MPI_PROC_NULL;

  ntimestep = 0;
  time = 0.0;
  stoptime = 0.0;
  ndump = 0;
  swap_freq = 0;

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
  pop_num_global = 0;
  num_variables = 0;

  best_fitness_global = 1.0e12;
  best_tree_global = NULL;
  best_pop_global = 0;

  thi = 0.0;
  tlo = 0.0;
  swap_freq = 0;

  temp_init = "linear";
}
/* ---------------------------------------------------------------------- */
AppGPPT::~AppGPPT()
{
  delete random;
  delete tree;
  delete fit;
  delete [] fitness;
  delete [] pop;
  delete [] var;
  for (int b = 0; b < buf_size; b++) delete [] buffer[b];
  delete [] buffer;
  delete [] tight_buf;
  for(int p = 0; p < pop_num; p++) delete [] tight_buf_pop[p];
  delete [] tight_buf_pop;
}
/* ---------------------------------------------------------------------- */
void AppGPPT::init()
{
  Node *current_root;
  double current_t;
  int pp;

  if(num_variables > 0 & range > 0)
    tree->init(const_lo, const_hi, num_variables, max_tree_depth);
  else {
    error->all("Number of variables or constant range not set.");
    return;
  }

  if(num_variables == fit->get_num_var()) fit->init();
  else
    error->all("Number of variables inconsistent in tree and fitness.");

  if(swap_freq == 0) error->all("Swap frequency is not set.");
  if (pop_num_global == 0) error->all("Number of populations not set.");
  if (pop_size == 0) error->all("Population size not set.");

  // distribute populations to procs.
  // offsets are calculated using integer math which 
  // automatically accounts for remainder populations.
  // local population counts are computed from offsets
  pop_offset = (me*pop_num_global)/nprocs;
  pop_num = ((me+1)*pop_num_global)/nprocs - pop_offset;

  if (pop_num <= 0) {
    error->all("Must have at least one population for each process");
  }

  pop = new Population[pop_num];
  current_t = thi;
  // Loop over global pops to preserve invariance to processor count 
  for(int p = 0; p < pop_num_global; p++){
    pp = pop_num_global-p-1 - pop_offset;
    if (pp >= 0 && pp < pop_num) {
      pop[pp].init(seed+p, pop_size, 
		   tree, fit, max_tree_depth, current_t);
    }
    if (strcmp(temp_init.c_str(),"log")==0) current_t = (tlo + current_t)/2; 
    else if (strcmp(temp_init.c_str(),"linear")==0)
      current_t -= (thi-tlo)/static_cast<double>(pop_num_global);
    else error->all("Illegal temperature distribution");
  }

  buf_size = 2;
  word_size = 16;
  for (int i = 0; i < max_tree_depth; i++) buf_size *= 2;
  buffer = new char*[buf_size];
  for (int b = 0; b < buf_size; b++) buffer[b] = new char[word_size];
  tight_buf = new char[buf_size*word_size];
  tight_buf_pop = new char*[pop_num];
  for(int p = 0; p < pop_num; p++){
    tight_buf_pop[p] = new char[buf_size*word_size];
  }

  if (me == 0) {
    if (screen) {
      fprintf
	(screen,"Step Time BestPop BestFitness");
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf
	(logfile,"Step Time BestPop BestFitness");
      fprintf(logfile,"\n");
    }
  }

  stats();
  
  // setup future calls to stats()
  
  stats_time = time + stats_delta;
  if (stats_delta == 0.0) stats_time = stoptime;
}
/* ---------------------------------------------------------------------- */
void AppGPPT::write_tex_header(Node *tree)
{
  if (me == 0) {
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
}
/* ---------------------------------------------------------------------- */
void AppGPPT::write_tree(int t)
{
  fprintf(screen,"######################## Tree %d############.\n",t);
  trees[t]->write(screen);
  fprintf(screen,"\n");
}

/* ---------------------------------------------------------------------- */
void AppGPPT::input(char *command, int narg, char **arg)
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
void AppGPPT::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  stoptime = time + atof(arg[0]);

  // init classes used by this app
  
  init();
  timer->init();

  // perform the run

  iterate();

  // final statistics

  Finish finish(spk);
}
/* ----------------------------------------------------------------------
   iterate on solver
------------------------------------------------------------------------- */
void AppGPPT::iterate()
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
    for(p = 0; p < pop_num; p++) pop[p].build_new_population();

    timer->stamp(TIME_SOLVE);

    //update mutation weights
    for(p = 0; p < pop_num; p++) pop[p].update_weights();
    //swap trees
    if (swap_freq > 0 && ntimestep%swap_freq == 0) {
      swap_trees(par);
      par = par?0:1;
    }

    if (ndump > 0 && ntimestep%ndump == 0) dump(); 
    
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
void AppGPPT::stats()
{
  int i;

  get_best_fitness();

  if (me == 0) {
    if (screen) {
      fprintf(screen,"%d %g %d %g ",ntimestep,time, best_pop_global, best_fitness_global);
      fprintf(screen,"\n \n");
    }
    if (logfile) {
      fprintf(logfile,"%d %g %d %g ",ntimestep,time, best_pop_global, best_fitness_global);
      fprintf(logfile,"\n \n");
    }

    if (screen){
      best_tree_global->write(screen);
      fprintf(screen, "\n \n");
      best_tree_global->write_tex(screen);
      fprintf(screen, "\n \n");
    }
    if (logfile){
      best_tree_global->write(logfile);
      fprintf(logfile, "\n \n");
      best_tree_global->write_tex(logfile);
      fprintf(logfile, "\n \n");
    }
  }

}
/* ----------------------------------------------------------------------
   proc 0 writes dump header
------------------------------------------------------------------------- */

void AppGPPT::dump_header()
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

// *** THIS DOES NOT WORK ANY MORE, EXCEPT MAYBE ON ONE PROCESSOR ***//
void AppGPPT::dump()
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
void AppGPPT::set_tree_type(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal tree_type command");
  max_tree_depth = atoi(arg[0]);

  const_lo = atof(arg[1]);
  const_hi = atof(arg[2]);
  range = const_hi - const_lo + 1;
  num_variables = atoi(arg[3]);
  if (me == 0) {
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
}
/* ---------------------------------------------------------------------- */
void AppGPPT::set_population(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal population command");
  pop_num_global = atoi(arg[0]);
  pop_size = atoi(arg[1]);
  if (me == 0) {
    if(screen){
      fprintf(screen,"Population number %d\n", pop_num_global);
      fprintf(screen,"Population size %d\n", pop_size);
    }
    if(logfile){
      fprintf(logfile,"Population number %d\n", pop_num_global);
      fprintf(logfile,"Population size %d\n", pop_size);
    }
  }
}
/* ---------------------------------------------------------------------- */
void AppGPPT::set_tempering(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal tempering command");

  temp_init.assign(arg[0]);

  if(!(temp_init == "log") & !(temp_init == "linear"))
    error->all("Illegal temperature distribution");

  tlo = atof(arg[1]);
  thi = atof(arg[2]);
  swap_freq = atoi(arg[3]);

  if (me == 0) {
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
}
/* ---------------------------------------------------------------------- */
void AppGPPT::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);

  if (me == 0) {
    if(screen){
      fprintf(screen,"Stats frequency = %g\n", stats_delta);
    }
    if(logfile){
      fprintf(logfile,"Stats frequency = %g\n", stats_delta);
    }
  }
}
/* ---------------------------------------------------------------------- */

// *** dump() DOES NOT WORK ANY MORE, EXCEPT MAYBE ON ONE PROCESSOR ***//
void AppGPPT::set_dump(int narg, char **arg)
{
  error->all("Dump command temporarily disabled");
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
void AppGPPT::set_fitness(int narg, char **arg)
{
  fit->set_fitness(narg, arg, me);
}
/* ---------------------------------------------------------------------- */
void AppGPPT::get_best_fitness()
{
  double current_best_fitness;
  char* inc_buf;
  int best_proc, best_pop;
  MPI_Status status;
  double best_fitness,best_fitness_list[nprocs];
  Node* best_tree;

  best_fitness = 1.0e12;
  best_pop = 0;

  for(int p = 0; p < pop_num; p++) {
    current_best_fitness = pop[p].get_best_fitness();
    if(current_best_fitness < best_fitness){
      best_fitness = current_best_fitness;
      best_pop = p;
    }
  }

  MPI_Allgather(&best_fitness,1,MPI_DOUBLE,&best_fitness_list,1,MPI_DOUBLE,world);

  best_fitness_global = 1.0e12;
  best_proc = 0;

  for(int iproc = 0; iproc < nprocs; iproc++) {
    if(best_fitness_list[iproc] < best_fitness_global){
      best_fitness_global = best_fitness_list[iproc];
      best_proc = iproc;
    }
  }

  if (best_proc != 0) {
    if (best_proc == me) {
      // Pack best into buffer
      best_tree = pop[best_pop].get_best_tree();
      pack_buffer(best_tree);
      // Send buffer to root process
      MPI_Send(tight_buf,buf_size*word_size,MPI_CHAR,0,0,world);
      // Send best_pop to root process
      best_pop += pop_offset;
      MPI_Send(&best_pop,1,MPI_INT,0,1,world);
    }

    if (me == 0) {
      MPI_Recv(tight_buf,buf_size*word_size,MPI_CHAR,best_proc,0,world,&status);
      best_tree_global = unpack_buffer();
      MPI_Recv(&best_pop_global,1,MPI_INT,best_proc,1,world,&status);
    }
  } else {
    best_tree_global = pop[best_pop].get_best_tree();
    best_pop_global = best_pop;
  }
}

/* ---------------------------------------------------------------------- */
void AppGPPT::swap_trees(int par)
{
  Population *curr_pop;
  double temp_send,temp_recv,temp_nbr;
  Node* test_tree;
  int send_nbr,recv_nbr;
  MPI_Status status;

  //populations pack buffers with best trees
  //neigbors receive
  for(int p = 0; p < pop_num; p++) {
    curr_pop = &pop[p];
    // pack best tree into tight_buf
    pack_buffer(curr_pop->get_best_tree());
    //determine send neighbor population index
    send_nbr = (par == 0)?p+1:p-1;

    if (send_nbr == -1) {
      MPI_Send(tight_buf,buf_size*word_size,MPI_CHAR,procdown,1,world);
      MPI_Recv(tight_buf_pop[pop_num-1],buf_size*word_size,MPI_CHAR,procup,1,world,&status);
      temp_send = curr_pop->get_temperature();
      MPI_Send(&temp_send,1,MPI_DOUBLE,procdown,2,world);
      MPI_Recv(&temp_recv,1,MPI_DOUBLE,procup,2,world,&status);
    } else if (send_nbr == pop_num) {
      MPI_Send(tight_buf,buf_size*word_size,MPI_CHAR,procup,3,world);
      MPI_Recv(tight_buf_pop[0],buf_size*word_size,MPI_CHAR,procdown,3,world,&status);
      temp_send = curr_pop->get_temperature();
      MPI_Send(&temp_send,1,MPI_DOUBLE,procup,4,world);
      MPI_Recv(&temp_recv,1,MPI_DOUBLE,procdown,4,world,&status);
    } else {
      memcpy(tight_buf_pop[send_nbr],tight_buf,buf_size*word_size);
    }
  }

  for(int p = 0; p < pop_num; p++) {
    curr_pop = &pop[p];
    // Skip over populations on the ends
    if (par == 1 && p+pop_offset == pop_num_global-1) continue;
    if (par == 0 && p+pop_offset == 0) continue;
    memcpy(tight_buf,tight_buf_pop[p],buf_size*word_size);
    //create tree from buffer and check acceptance
    test_tree = unpack_buffer();
    //determine recv neighbor population index
    recv_nbr = (par == 1)?p+1:p-1;
    if (recv_nbr == -1) {
      temp_nbr = temp_recv;
    } else if (recv_nbr == pop_num) {
      temp_nbr = temp_recv;
    } else {
      temp_nbr = pop[recv_nbr].get_temperature();
    }
    int swapflag;
    double testfit,bestfit;
    curr_pop->swap_nbr(par, temp_nbr, test_tree);
  }
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
void AppGPPT::pack_buffer(Node *root)
{
  int buf_cnt;
  int n;

  // Initialize tight buffer to \x1
  memset(tight_buf,'\x1',buf_size*word_size);
  //   printf("\npack buffer before %d \n",n);
  //   write_tight_buf();

  // First pack tree into buffer
  // Initialize each word in buffer to \x0
  for(int w = 0; w < buf_size; w++) memset(buffer[w],'\x0',word_size);

  // Copy tree into buffer using recursion on root node
  n = 0;
  root->buffer(buffer, n);

  // Then pack buffer into tight buffer
  buf_cnt = 0;
  for(int w = 0; w < n; w++) {
    for (int i = 0; i < word_size; i++){
      tight_buf[buf_cnt] = buffer[w][i];
      buf_cnt ++;
    }
  }
  tight_buf[buf_cnt] = endbufchar;

//   printf("\npack buffer after %d \n",n);
//   write_tight_buf();
}

/* ---------------------------------------------------------------------- */
Node* AppGPPT::unpack_buffer()
{
  int buf_cnt;
  char current_ch;
  int w = 0;
  int n = 0;

//   printf("\nunpack buffer before %d \n",n);
//   write_tight_buf();

  buf_cnt = 0;
  while(1==1){
    for (int i = 0; i < word_size; i++){
      current_ch = tight_buf[buf_cnt];
      buffer[w][i] = current_ch;
      buf_cnt ++;
      if(current_ch == endbufchar){
	return tree->from_buffer(buffer, n);
      }
    }
    w++;
  }

}

void AppGPPT::write_tight_buf() {
  int buf_cnt = 0;
  for (int i=0;i<buf_size;i++) {
    printf("\n word %d ",i);
    for (int j=0;j<word_size;j++) {
      printf("%x ",tight_buf[buf_cnt++]);
    }
  }
  cout << endl;
}
