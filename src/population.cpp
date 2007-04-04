/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
#include <iostream>
#include <sstream>
#include "random_park.h"
#include "node.h"
#include "tree.h"
#include "population.h"
#include <cmath>
#include "string.h"
#include "plus_node.h"
#include "minus_node.h"
#include "times_node.h"
#include "divide_node.h"
#include "const_node.h"
#include "var_node.h"
#include "error.h"

/* ---------------------------------------------------------------------- */
using namespace std;
using namespace SPPARKS;

Population::Population() : SysPtr(spk)
{
}
/* ---------------------------------------------------------------------- */
Population::~Population()
{
  for(int t=0;t<ntrees;t++) roots[t]->clear();
  delete [] roots;
  delete [] fitness;
  for (int b = 0; b < buf_size; b++) delete [] buffer[b];
  delete [] buffer;
  delete [] tight_buf;
}
/* ---------------------------------------------------------------------- */
void Population::init(int seed, int num_in, Tree * tree_in,
		      Fitness *fit_in, int max_depth, double temp_in)
{
  if (random != NULL)delete random;
  random = new RandomPark(seed);

  tree = tree_in;
  ntrees = num_in;
  temperature = temp_in; beta = 1.0/temp_in;
  fit = fit_in;

  fitness = new double[ntrees];
  roots = new Node*[ntrees];

  best_fit = 1e12;
  best_index = 0;

  for(int t=0;t<ntrees;t++){
    roots[t] = tree->build_tree(max_depth);
    fitness[t] = fit->compute(roots[t]);
    if (fitness[t] < best_fit) {
      best_fit = fitness[t];
      best_index = t;
    }
  }

  //allocate buffers
  buf_size = 2;
  word_size = 16;
  for (int i = 0; i < max_depth; i++) buf_size *= 2;
  buffer = new char*[buf_size];
  for (int b = 0; b < buf_size; b++) buffer[b] = new char[word_size];
  tight_buf = new char[buf_size*word_size];

  set_weights();

  swap_up = 0;
  swap_down = 0;
}

/* ---------------------------------------------------------------------- */
void Population::build_new_population()
{
  Node *tree1;
  Node *tree2;
  int mut_oper = 0;
  int partner;
  int t = 0;
  double current_fit = 0;

  //create new trees from old trees
  while (t < ntrees){
    //select operator
    accepted[MUT_SWAP+1]++;
    mut_oper = linear_select_mutation();
    //next tree
    tree1 = tree->copy(roots[t]); t++;
    //build new tree(s)
    if      (mut_oper == MUT_CONST)
      {tree->mutate_constant(tree1, temperature);}
    else if (mut_oper == MUT_VAR){tree->mutate_variable(tree1);}
    else if (mut_oper == MUT_OPER){tree->mutate_operator(tree1);}
    else if (mut_oper == MUT_BRANCH){tree->mutate_branch(tree1);}
    else if ((mut_oper == MUT_CROSS ||mut_oper == MUT_SWAP) & t < ntrees - 1){
      tree2 = tree->copy(roots[t]); 
      t++;
      if(mut_oper == MUT_CROSS) 
	tree->crossover(tree1, tree2, tree->max_depth);
      if (accept(roots[t-2],tree2, current_fit)){
	accepted[mut_oper]++;
	roots[t-2]->clear();
	roots[t-2] = tree2;
	//update best fit
	if (current_fit < best_fit) {
	  best_fit = current_fit;
	  best_index = t-2;
	}
	//update fitness
	fitness[t-2] = current_fit;
      }
      else tree2->clear();
    }

    //check acceptance and replace if needed
    if(accept(roots[t-1],tree1, current_fit)){
      accepted[mut_oper]++;
      roots[t-1]->clear();
      roots[t-1] = tree1;
      fitness[t-1] = current_fit;
      if (current_fit < best_fit){
	best_fit = current_fit;
	best_index = t-1;
      }
    }
    else tree1->clear();
  }
}
/* ---------------------------------------------------------------------- */
double Population::get_stats(double &sigma)
{
  double mean = 0.0;
  sigma = 0.0;
  double cur_fit;

  for (int t = 0; t < ntrees; t++){
    cur_fit = fitness[t];
    mean += cur_fit;
    sigma += cur_fit*cur_fit;
  }
  sigma /= static_cast<double>(ntrees);
  mean /= static_cast<double>(ntrees);

  sigma -= mean*mean;
  sigma = sqrt(sigma);

  return mean;
}
/* ---------------------------------------------------------------------- */
void Population::set_weights()
{
  double sum =0.0;
  int i;
  //set weights

  for (i = 0; i < MUT_SWAP+1; i++) weight[i] = 1.0;
  for (i = 0; i < MUT_SWAP+1; i++) sum += weight[i];
  for (i = 0; i < MUT_SWAP+1; i++) weight[i] /= sum;

  for(i = 0; i < MUT_SWAP+2; i++)
    accepted[i]=0;
}
/* ---------------------------------------------------------------------- */
void Population::update_weights()
{
  double acc_ratio[MUT_SWAP+1];
  int i;
  double offset;

  for (i = 0; i < MUT_SWAP+1; i++) {
    acc_ratio[i] = static_cast<double>(accepted[i])/
      static_cast<double>(accepted[MUT_SWAP+1]);
    offset = acc_ratio[i] - 0.5*weight[i];
    weight[i] += offset;
    if(weight[i] <= 0.1) weight[i] = 0.2;
  }

  double sum = 0.0;

  for (i = 0; i < MUT_SWAP+1; i++) sum += weight[i];
  for (i = 0; i < MUT_SWAP+1; i++) weight[i] /= sum;
  for (i = 0; i < MUT_SWAP+2; i++) accepted[i]=0;
}
/* ---------------------------------------------------------------------- */
void Population::write_best_tree(FILE * dest)
{
  fprintf(dest, "Best tree: \n");
  roots[best_index]->write(dest);
  fprintf(dest, "Fit: %g \n", fitness[best_index]);
  fprintf(dest, "Tree value: %g \n",fit->tree_value(roots[best_index]));
}
/* ---------------------------------------------------------------------- */
Node *Population::get_best_tree()
{
  return roots[best_index]; 
}
/* ---------------------------------------------------------------------- */
Node *Population::get_tree(int t)
{
  return roots[t]; 
}
/* ---------------------------------------------------------------------- */
bool Population::accept(Node *tree1, Node *tree2, double &best_fit)
{
  double fit1, fit2;

  fit1 = fit->compute(tree1);
  fit2 = fit->compute(tree2);

  if(beta*(fit2-fit1) < log(random->uniform())) {
    best_fit = fit2;
    return true;
  }
  best_fit = fit1;
  return false;
}
/* ----------------------------------------------------------------------
   Sample mutation distribution with partial sums
   ------------------------------------------------------------------------- */
int Population::linear_select_mutation()
{
  int g = 0;
  double compare = random->uniform();
  double partial = 0.0;

  while (1==1){
    partial += weight[g];
    if (partial > compare) return g;
    g++;
  }
  error->all("Mutation operator selection failed.");

}
/* ---------------------------------------------------------------------- */
void Population::pack_buffer(Node *root)
{
  int n = 0;
  char endc = '#';

  buf_cnt = 0;
  root->buffer(buffer, n);

  for(int w = 0; w < n; w++)
    for (int i = 0; i < word_size; i++){
      tight_buf[buf_cnt] = buffer[w][i];
      buf_cnt ++;
    }
  tight_buf[buf_cnt] = endc;
}
/* ---------------------------------------------------------------------- */
void Population::unpack_buffer(char *buf)
{
  char endc = '#';
  char current_ch;
  buf_cnt = 0;
  int w = 0;

  while(1==1){
    for (int i = 0; i < word_size; i++){
      current_ch = buf[buf_cnt];
      buffer[w][i] = current_ch;
      buf_cnt ++;
      if(current_ch == endc){
	return;
      }
    }
    w++;
  }
}
/* ---------------------------------------------------------------------- */
void Population::swap_nbr(int par, double t_in)
{
  int n = 0;
  double beta_in = 1.0/t_in;

  Node *test_tree = tree->from_buffer(buffer, n);

  double fit1 = 0.0;

  fit1 = fit->compute(test_tree);

  if( (fit1-best_fit) * (beta_in-beta) < log(random->uniform())) {
    //    cout <<"best_fit = "<<best_fit;
    best_fit = fit1;
    //    cout << " new fit = "<<best_fit<<" par "<<par<<endl;
    roots[best_index]->clear();
    roots[best_index] =  test_tree;
    if (par == 0) swap_up++;
    else if(par == 1) swap_down ++;
  }
}
