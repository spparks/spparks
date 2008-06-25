/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "state.h"
#include "random_park.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

State::State(int size_in, int seed) 
{
  size = size_in;
  nvars_current = 0;

  rnd = new RandomPark(12345);
}

/* ---------------------------------------------------------------------- */

State::~State()
{
}

/* ---------------------------------------------------------------------- */
void State::init_var(int narg, char **arg)
{

  double lo, hi;
  double cd = 0.0;
  int ci = 0;
  int v;

  double *dp;
  int *ip;
  int n;

  //  strcpy(dist_type,tokens[0]);

  n = get_attribute(arg[0]);


  if(data_type[n] == DOUBLE_TYPE) dp = (double *) attr[n];
  else if(data_type[n] == INT_TYPE) ip = (int *) attr[n];
  // else cout <<"variable "<<n<<" type not set."<<endl;

  //  cout << "distribution type = "<< dist_type<<endl;

  if(strcmp(arg[1],"constant")==0){

    //  cout <<"constant = "<<c<<endl;
    if(data_type[n] == DOUBLE_TYPE){
      cd = atof(arg[2]);
      for(v = 0; v < size; v++) dp[v] = cd;
    }
    else if(data_type[n] == INT_TYPE) {
      ci = atoi(arg[2]);
      for(v = 0; v < size; v++) ip[v] = ci;
    }
  }
  else if(strcmp(arg[1],"linear")==0){
    lo =  atof(arg[2]);
    hi =  atof(arg[3]);
  }
  else if(strcmp(arg[1],"random")==0){
    lo =  atof(arg[2]);
    hi =  atof(arg[3]);
  }
  //transfer check to the app
  //  else error->one("Unknown distribution type in init.");

}

/* ---------------------------------------------------------------------- */

void State::add_attribute(char *name_in, int data_type_in)
{
  bool bl;
  int i;
  double db;

//   cout <<"setting variable"<<endl;
//   cout <<"data type = "<< data_type_in<<endl;
//   cout <<"name = "<<name_in<<endl;
//   cout<< "index in the array = "<< nvars_current<<endl;

  data_type[nvars_current] = data_type_in;
  if(data_type_in == BOOL_TYPE) step[nvars_current] = sizeof(bool);
  if(data_type_in == INT_TYPE) step[nvars_current] = sizeof(int);
  if(data_type_in == DOUBLE_TYPE) step[nvars_current] = sizeof(double);

  int nsize = size*step[nvars_current];
  void* temp = (void*)malloc(nsize);
  strcpy(name[nvars_current],name_in);

  attr[nvars_current] = temp;
  nvars_current ++;
}

/* ---------------------------------------------------------------------- */

void *State::get_attribute(char *name_in, int &type_out)
{
  int cnt = 0;

  while (strcmp(name[cnt],name_in)!=0 & cnt < NVARS_MAX) cnt++;
  if(cnt<NVARS_MAX) {
    type_out = data_type[cnt];
    return attr[cnt];
  }
  //transfer check to app
  //  else cout << "failed to find variable "<<name_in<<endl;

  return NULL;
}

/* ---------------------------------------------------------------------- */

int State::get_attribute(char *name_in)
{
  int cnt = 0;

  while (strcmp(name[cnt],name_in)!=0 & cnt < NVARS_MAX) cnt++;
  if(cnt<NVARS_MAX) {
    return cnt;
  }
  //transfer check to app
  //  else error->one("Illegal state command");

  return -1;
}
