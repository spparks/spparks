/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

/*----------------------------------------------------------------------

----------------------------------------------------------------------*/

#ifndef APP_GPPT_H
#define APP_GPPT_H

#include "app.h"
#include "node.h"
#include <string>
#include "solve_gppt.h"
#include "fitness.h"
#include "population.h"

using namespace std;
namespace SPPARKS {

class AppGppt : public App {
 public:
  AppGppt(class SPK *, int, char **);
  ~AppGppt();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

  Node **trees;

  //SolveGPPT *solve;

 private:

  class RandomPark *random;

  Tree *tree;
  Fitness *fit;
  Population *pop;

  //parallel
  int me,nprocs;
  int procup, procdown;    // neighbor procs
  int pop_per_proc;

  //time keeping
  int ntimestep;
  double time,stoptime;


  //stats timing
  double stats_time, stats_delta;

  void iterate();
  void stats();

  //commands
  void set_stats(int, char **);
  void set_tree_type(int, char **);
  void set_population(int, char **);
  void set_tempering(int, char **);
  void set_fitness(int, char **);
  void set_dump(int, char **);

  //output
  void write_tex_header(Node *);
  void write_tree(int);
  void dump_header();
  void dump();

  FILE *fp;
  int ndump;
  int dpop1, dpop2;

  //tree parameters
  int max_tree_depth;
  double const_lo;
  double const_hi;
  double range;
  int num_variables;


  //population parameters
  int pop_size;
  int pop_num;

  //tempering parameters
  double temperature;
  string temp_init;
  double thi, tlo;
  int swap_freq;
  void swap_trees(int);

  //tree buffering
  char** buffer;
  int buf_size;
  int word_size;
  char *tight_buf;

  //fitness variables
  double *fitness;
  double *var;

  double best_fitness;
  Node *best_tree;
  int best_pop;
  double get_best_fitness();

  //random seed
  int seed;

};

}

#endif
