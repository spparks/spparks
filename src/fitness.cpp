/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
#include <iostream>
#include "random_park.h"
#include "node.h"
#include "fitness.h"
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
//#include "spk.h"

/* ---------------------------------------------------------------------- */
using namespace std;
using namespace SPPARKS;

Fitness::Fitness(SPK *spk) : SysPtr(spk)
{
  nvar = 0;

  random = new RandomPark(321123);

}
/* ---------------------------------------------------------------------- */
Fitness::~Fitness()
{
  delete [] var;
  for(int n = 0; n < nsnaps; n++) delete [] snaps[n];
  delete [] snaps;
  delete [] energy;
  delete [] test_e;
  delete [] npairs;
  delete [] hits;
  delete [] weight;
}
/* ---------------------------------------------------------------------- */
void Fitness::init()
{
  if(nvar > 0){
    var = new double[nvar];
  }
  else error->all("Number of variables not set.");

}
/* ---------------------------------------------------------------------- */
double Fitness::compute(Node * root)
{
  int n;
  int p;
  double least_sq = 0.0;
  double diff;
  double test;
  double max = -1.0e33;
  int imax = 0;

  for(n=0;n<nsnaps;n++){
    test_e[n] = 0;
    for(p=0;p<npairs[n];p++){
      var[0] = snaps[n][p];
      test_e[n] += root->go(var);
    }  
  }

  for(n=0;n<nsnaps;n++) {
    diff = energy[n]-test_e[n];
    test = abs(diff);
    if (test > max) {max = test; imax = n;}
    least_sq += weight[n]*diff*diff;
  }

  hits[imax]++;
  hits[nsnaps]++;
  if(hits[nsnaps] == 100000) {
//     cout <<"energy: ";
//     for(n=0;n<nsnaps;n++)cout <<energy[n]<<" ";
//     cout <<endl;
//     cout <<"test_e: ";
//     for(n=0;n<nsnaps;n++)cout <<test_e[n]<<" ";
//     cout <<endl;
//     cout <<"diff: ";
//     for(n=0;n<nsnaps;n++)cout <<abs(energy[n]-test_e[n])<<" ";
//     cout <<endl;
//     cout <<"hits: ";
//     for(n=0;n<nsnaps+1;n++)cout <<hits[n]<<" ";
//     cout <<endl;

    update_weight(100000);
//     cout <<"weights: ";
//     for(n=0;n<nsnaps;n++)cout <<weight[n]<<" ";
//     cout <<endl;

  }

  return sqrt(least_sq);
}
/* ---------------------------------------------------------------------- */
void Fitness::update_weight(int sum)
{
  int n;
  double sumd = static_cast<double>(sum);
  
  for(n=0;n<nsnaps;n++)
    if(hits[n] < 1000) 
      {sumd += 10000 - hits[n];hits[n] = 10000;}
  for(n=0;n<nsnaps;n++)
    weight[n] = static_cast<double>(hits[n])/sumd;
  for(n=0;n<nsnaps+1;n++) hits[n] = 0;

}
/* ---------------------------------------------------------------------- */
void Fitness::set_fitness(int narg, char **arg, int me)
{
  if (narg != 2) error->all("Illegal fitness command");

  nvar = atoi(arg[0]);
  fitfile_name.assign(arg[1]);

  if(me == 0){
    fitfile = fopen(fitfile_name.c_str(),"r");
    if (fitfile == NULL) 
      error->all("Cannot open fitness data file");
    if(screen){
      fprintf(screen,"Sensor variables in fitness = %d \n", nvar);
      fprintf(screen,"Reading fitness data from ");
      fprintf(screen, fitfile_name.c_str());
      fprintf(screen,"\n");
    }
    if(logfile){
      fprintf(logfile,"Sensor variables in fitness = %d\n", nvar);
      fprintf(logfile,"Reading fitness data from ");
      fprintf(logfile, fitfile_name.c_str());
      fprintf(logfile,"\n");
    }
  }

  //  read_fit_file(me);
  //**********************************************
  //this is a temporary hack manufacturing fitness data
  //to replace, uncomment the line above

  nsnaps = 25;
  
  snaps  = new double *[nsnaps];
  energy = new double[nsnaps];
  test_e = new double[nsnaps];
  npairs = new int[nsnaps];
  hits   = new int[nsnaps+1];
  weight = new double[nsnaps];

  for(int n=0;n<nsnaps+1;n++) hits[n] = 0;
  for(int n=0;n<nsnaps;n++) weight[n]=1.0/nsnaps;

  for (int n = 0; n < nsnaps; n++){
    npairs[n] = 6;
    snaps[n] = new double[npairs[n]];
  }

  for (int n = 0; n < nsnaps; n++)
    for (int p = 0; p < npairs[n]; p++){
      snaps[n][p] = get_distance();
    }
  //end of temporary hack
  //****************************************
  
  for (int n = 0; n < nsnaps; n++){
    energy[n] = 0.0;
    for (int p = 0; p < npairs[n]; p++){
      double over_d_six = pow(1.0/snaps[n][p],6); 
      energy[n] += (over_d_six-1.0)*over_d_six;
    }
  }
}
  /* ---------------------------------------------------------------------- */
double Fitness::get_distance()
{
  double lo = 0.45;
  double hi = 4.0;
  double range = 3.55;
  double vol = (hi*hi*hi-lo*lo*lo)/3.0;
  double shell = vol;
  double dist;

  while (shell/vol > random->uniform()){
    dist = lo + range*random->uniform();
    shell = dist*dist;
  }

  return dist;
}
  /* ---------------------------------------------------------------------- */
void Fitness::read_fit_file(int me)
{
  int cnt; 

  if(me == 0) fscanf(fitfile,"NSNAPS %d \n", &nsnaps);

  MPI_Bcast(&nsnaps,1,MPI_INT,0,world);

  snaps = new double *[nsnaps];
  energy = new double[nsnaps];
  test_e = new double[nsnaps];
  npairs = new int[nsnaps];
  hits   = new int[nsnaps+1];
  weight = new double[nsnaps];

  for(int n=0;n<nsnaps+1;n++) hits[n] = 0;
  for(int n=0;n<nsnaps;n++) weight[n]=1.0/nsnaps;


  for (int n = 0; n < nsnaps; n++){
    if(me == 0) {
      fscanf(fitfile,"SNAP %d ENERGY %lg PAIRS %d \n", 
	     &cnt, &energy[n], &npairs[n]);
    }

//     if(me == 0) fprintf(screen, "SNAP %d ENERGY %lg PAIRS %d \n",
// 	    cnt, energy[n], npairs[n]);

    MPI_Bcast(&energy[n],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&npairs,1,MPI_INT,0,world);

    snaps[n] = new double[npairs[n]];

    if(me == 0)
      for (int p = 0; p < npairs[n]; p++){
	fscanf(fitfile,"%lg \n", &snaps[n][p]);
	//	fprintf(screen,"%lg \n",snaps[n][p]);
      }

    MPI_Bcast(&snaps[n],npairs[n],MPI_DOUBLE,0,world);
  }


  if(me == 0) {
    if (nsnaps == 0) 
      error->all("No training set read from file.");
    if(screen){
      fprintf(screen,"Number of snapshots in training set = %d \n", nsnaps);
    }
    if(logfile){
      fprintf(logfile,"Number of snapshots in training set = %d \n", nsnaps);
    }
  }
}
