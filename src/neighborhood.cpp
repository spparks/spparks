#include "neighborhood.h"
#include <iostream>
#include "random_park.h"
#include "state.h"


using namespace std;
using namespace SPPARKS;

/* ---------------------------------------------------------------------- */
Neighborhood::Neighborhood(char *name_in, State *state_in)
{
  nbname.assign(name_in);
  state = state_in;
  size = state->get_size();
}
/* ---------------------------------------------------------------------- */
Neighborhood::~Neighborhood()
{

}
/* ---------------------------------------------------------------------- */
void Neighborhood::init(char *type)
{
  string input(type);


}
/* ---------------------------------------------------------------------- */
