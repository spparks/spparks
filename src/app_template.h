/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_TEMPLATE_H
#define APP_TEMPLATE_H

#include "app.h"
#include "node.h"
#include "state.h"
#include "expression.h"
#include <vector>

using namespace std;
namespace SPPARKS {
  
  class AppTemplate : public App {
  public:
    AppTemplate(class SPK *,int, char **);
    virtual ~AppTemplate();
    void init();
    void input(char *, int, char **);
    void run(int, char **);
    void stats(char *);
    void stats_header(char *);

  private:
    int ntimestep;
    double time,stoptime;
    double *propensity;        // propensity of each event
    int nevents;
    
    void iterate();
    
    void set_stats(int, char **);
    void set_state(int, char **);
    void set_variable(int, char **);
    void init_variable(int, char **);
    void set_expression(int, char **);
    
    State *st;
    Tree *tr;
    Xpression **xpr;
    int nxpr;
    vector<Xpression *> expr;
  };
  
}

#endif
