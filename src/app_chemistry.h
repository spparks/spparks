/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_APP_CHEMISTRY_H
#define SPK_APP_CHEMISTRY_H

#include "app.h"

namespace SPPARKS_NS {

class AppChemistry : public App {
 public:
  AppChemistry(class SPPARKS *, int, char **);
  virtual ~AppChemistry();
  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();

 private:
  double volume;                 // Gillespie volume

  int nevents;                   // # of reactions performed
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

  void stats(char *);
  void stats_header(char *);

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
