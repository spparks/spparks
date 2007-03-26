/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_SURF_H
#define APP_SURF_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppSurf : public App {
 public:
  AppSurf(class SPK *, int, char **);
  ~AppSurf();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

 protected:
  int me,nprocs;
  int nlattice,seed;
  double strain;
  double temperature;
  int nstats,ndump;
  double rate_deposit;
  double cutoff;
  int nsteps;

  int ntimestep;
  FILE *fp;
  double hop_distance;
  double attempt_frequency;
  int stats_next,dump_next;
  double xlo,xhi,xprd;
  double zlo,zhi;
  double time;
  double energy;
  double cutsq;
  double skin,cutneigh,cutneighsq;
  double cutcount,cutcountsq;
  double epsilon,sigma;
  double lj1,lj2,lj3,lj4,offset;

  class RandomPark *random;

  int maxatomneigh,maxneigh;
  int **firstneigh;
  int *numneigh;
  int *neigh;

  struct OneAtom {
    double x,z;
    double fx,fz;
    int id,type;
  };

  int nlocal,nghost,maxatom;
  OneAtom *atoms;

  struct OneEvent {
    int iatom,style;
    double x,z;
  };

  int nevents,maxevent;
  OneEvent *events;
  double *rates;

  void iterate();
  int count_neigh(int);
  double zmax();
  double zrelax(int, double *, double *, double *);
  void add_atom(int, int, int, double, double);
  void add_event(int, int, double, double, double);
  double find_barrier(int, double, double *, double *);
  void relax(int);
  double etotal();
  double eng_force(int, int *);
  void ghost_comm();
  void pbc();
  void neighbor(int, int *);
  void neighbor_z(int, double);
  double dbrent(int, double, double, double, double, double *);

  void sd(int, int *);
  void cg(int, int *);
  int linemin_scan(int, int *, double *, double &, double, double,
		   double &, int &);
  int linemin_secant(int, int *, double *, double &, double, double,
		     double &, int &);

  void stats();
  void dump();

  void set_temperature(int, char **);
  void set_potential(int, char **);
  void set_rates(int, char **);
  void set_stats(int, char **);
  void set_dump(int, char **);
};

}

#endif
