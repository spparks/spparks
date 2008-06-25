/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "population.h"
#include "node.h"
#include "random_park.h"
#include "fitness.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Population::Population() : SysPtr(spk)
{
}

/* ---------------------------------------------------------------------- */

Population::~Population()
{
  for(int t=0;t<ntrees;t++) roots[t]->clear();
  delete [] roots;
  delete [] fitness;
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
    roots[t] = tree->build_tree(random,max_depth);
    fitness[t] = fit->compute(roots[t]);
  }

  calc_best_fit();
  set_weights();

  swap_up = 0;
  swap_down = 0;
}

/* ---------------------------------------------------------------------- */

void Population::calc_best_fit(){
  best_fit = 1e12;
  best_index = 0;
  for(int t=0;t<ntrees;t++)
    if (fitness[t] < best_fit) {
      best_fit = fitness[t];
      best_index = t;
    }
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
      {tree->mutate_constant(random,tree1, temperature);}
    else if (mut_oper == MUT_VAR){tree->mutate_variable(random,tree1);}
    else if (mut_oper == MUT_OPER){tree->mutate_operator(random,tree1);}
    else if (mut_oper == MUT_BRANCH){tree->mutate_branch(random,tree1);}
    else if ((mut_oper == MUT_CROSS ||mut_oper == MUT_SWAP) & t < ntrees - 1){
      tree2 = tree->copy(roots[t]); 
      t++;
      if(mut_oper == MUT_CROSS) 
	tree->crossover(random,tree1, tree2, tree->max_depth);
      if (accept(roots[t-2],tree2, current_fit)){
	accepted[mut_oper]++;
	roots[t-2]->clear();
	roots[t-2] = tree2;
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
    }
    else tree1->clear();
  }
  //update best fit
  calc_best_fit();
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

bool Population::accept(Node *tree1, Node *tree2, double &fit_in)
{
  double fit1, fit2;

  fit1 = fit->compute(tree1);
  fit2 = fit->compute(tree2);

  if(beta*(fit1-fit2) > log(random->uniform())) {
    fit_in = fit2;
    return true;
  }
  fit_in = fit1;
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

void Population::swap_nbr(int par, double t_in, Node *test_tree)
{
  int n = 0;
  double beta_in = 1.0/t_in;

//   calc_best_fit();

  double fit1 = 0.0;

  fit1 = fit->compute(test_tree);

  if( (fit1-best_fit) * (beta_in-beta) > log(random->uniform())) {
    best_fit = fit1;
    roots[best_index]->clear();
    roots[best_index] = test_tree;
    if (par == 0) swap_up++;
    else if(par == 1) swap_down ++;
  }
}
