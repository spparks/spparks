/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_CHEMISTRY_H
#define APP_CHEMISTRY_H

#include "app.h"

namespace SPPARKS {

class AppChemistry : public App {
 public:
  AppChemistry(class SPK *, int, char **);
  virtual ~AppChemistry();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);

 private:
  int ntimestep;
  double time,stoptime;

  double volume;                 // Gillespie volume

  int nspecies;                  // # of unique species
  char **sname;                  // ID of each species

  int nreactions;                // # of user defined reactions
  char **rname;                  // ID of each reaction
  int *nreactant;                // nreactant[I] = # of reactants of reaction I
  int **reactants;               // reactants[I][J] = particle species of Jth
				 //   reactant of reaction I
  int *nproduct;                 // nproduct[I] = # of products of reaction I
  int **products;                // products[I][J] = particle species of Jth
				 //   product of reaction I
  double *rate;                  // rate[I] = input rate for reaction I

  int *pcount;               // counts for each species
  int *ndepends;             // # of reactions that depend on each reaction
  int **depends;             // i,j = jth reaction that depends on ith reaction
  double *propensity;        // propensity of each reaction

  double factor_zero;        // conversion factor for different reactions
  double factor_dual;

  int *rcount;               // statistics
  double stats_time,stats_delta;

  void iterate();
  void stats();

  void set_count(int, char **);
  void add_species(int, char **);
  void add_reaction(int, char **);
  void set_volume(int, char **);
  void set_stats(int, char **);

  int find_reaction(char *);
  int find_species(char *);
  void build_dependency_graph();
  double compute_propensity(int);
};

}

#endif
