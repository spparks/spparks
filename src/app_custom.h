/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_CUSTOM_H
#define APP_CUSTOM_H

#include "app_lattice.h"
#include "expression.h"
#include "neighborhood.h"
#include "event.h"
#include <map>
#include <vector>

namespace SPPARKS_NS {

class AppCustom : public AppLattice {
 public:
  AppCustom(class SPPARKS *, int, char **);
  ~AppCustom();

  double site_energy(int);
  void site_pick_random(int, double);
  void site_pick_local(int, double);
  double site_propensity(int, int);
  void site_event(int, int);
  void site_clear_mask(char *, int);

 private:
  char **name_int;
  char **name_dbl;
  int n_int_var;
  int n_dbl_var;
  int *sites;

  void input_app(char *, int, char **);
  void set_variable(int, char **);
  int get_variable(char *, int&);
  void init_variable(int, char **);
  void set_expression(int, char **);
  void set_neighborhood(int, char **);
  void set_event(int, char **);
  void set_propensity(int, char **);
  void set_energy(int, char **);

  std::map<int,int> hash;
  std::map<int,int>::iterator loc;

  class Tree *tr;

  std::vector<Expression *> expr;
  class Expression *get_expression(char *);
  std::vector<Neighborhood *> nbhd;
  class Neighborhood *get_neighborhood(char *);
  std::vector<Event *> events;
  class Event *get_event(char *);

  int *update_numneigh;
  int **update_neighbor;

  class Expression *energyxp;

  };

}

#endif
