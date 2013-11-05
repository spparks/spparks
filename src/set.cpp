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

#include "spktype.h"
#include "stdlib.h"
#include "string.h"
#include "set.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "region.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{IARRAY,DARRAY,X,Y,Z,XYZ,ID};
enum{VALUE,RANGE,UNIQUE,DISPLACE};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,TAGINT};

/* ---------------------------------------------------------------------- */

Set::Set(SPPARKS *spk) : Pointers(spk) {}

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (app->sites_exist == 0) error->all(FLERR,"Set command before sites exist");

  if (narg < 2) error->all(FLERR,"Illegal set command");

  int lhs,rhs;

  if (strcmp(arg[0],"site") == 0) {
    lhs = IARRAY;
    siteindex = 0;
    if (app->iarray == NULL)
      error->all(FLERR,"Setting a quantity application does not support");
  } else if (arg[0][0] == 'i') {
    lhs = IARRAY;
    siteindex = atoi(&arg[0][1]);
    if (siteindex < 1 || siteindex > app->ninteger)
      error->all(FLERR,"Setting a quantity application does not support");
    siteindex--;
  } else if (arg[0][0] == 'd') {
    lhs = DARRAY;
    siteindex = atoi(&arg[0][1]);
    if (siteindex < 1 || siteindex > app->ndouble)
      error->all(FLERR,"Setting a quantity application does not support");
    siteindex--;
  } else if (strcmp(arg[0],"x") == 0) {
    lhs = X;
  } else if (strcmp(arg[0],"y") == 0) {
    lhs = Y;
  } else if (strcmp(arg[0],"z") == 0) {
    lhs = Z;
  } else if (strcmp(arg[0],"xyz") == 0) {
    lhs = XYZ;
  } else error->all(FLERR,"Illegal set command");

  int iarg;
  if (strcmp(arg[1],"value") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal set command");
    rhs = VALUE;
    if (lhs == IARRAY) ivalue = atoi(arg[2]);
    else if (lhs == DARRAY) dvalue = atof(arg[2]);
    else error->all(FLERR,"Illegal set command");
    iarg = 3;
  } else if (strcmp(arg[1],"range") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal set command");
    rhs = RANGE;
    if (lhs == IARRAY) {
      ivaluelo = atoi(arg[2]);
      ivaluehi = atoi(arg[3]);
    } else if (lhs == DARRAY) {
      dvaluelo = atof(arg[2]);
      dvaluehi = atof(arg[3]);
    } else error->all(FLERR,"Illegal set command");
    iarg = 4;
  } else if (strcmp(arg[1],"unique") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal set command");
    rhs = UNIQUE;
    iarg = 2;
  } else if (strcmp(arg[1],"displace") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal set command");
    rhs = DISPLACE;
    if (lhs != X && lhs != Y && lhs != Z && lhs != XYZ)
      error->all(FLERR,"Illegal set command");
    dvalue = atof(arg[2]);
    iarg = 3;
  } else error->all(FLERR,"Illegal set command");
    
  // parse optional args

  fraction = 1.0;
  regionflag = 0;
  loopflag = 1;
  ncondition = 0;
  cond = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"fraction") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      fraction = atof(arg[iarg+1]);
      if (fraction <= 0.0 || fraction > 1.0) 
        error->all(FLERR,"Illegal set command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion < 0) error->all(FLERR,"Set command region ID does not exist");
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strcmp(arg[iarg+1],"all") == 0) loopflag = 1;
      else if (strcmp(arg[iarg+1],"local") == 0) loopflag = 0;
      else error->all(FLERR,"Illegal set command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"if") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      cond = (Condition *) memory->srealloc(cond,
					    (ncondition+1)*sizeof(Condition),
					    "set:cond");
      if (strcmp(arg[iarg+1],"id") == 0) {
	cond[ncondition].lhs = ID;
	cond[ncondition].type = TAGINT;
	cond[ncondition].stride = 1;
      } else if (strcmp(arg[iarg+1],"x") == 0) {
	cond[ncondition].lhs = X;
	cond[ncondition].type = DOUBLE;
	cond[ncondition].stride = 3;
      } else if (strcmp(arg[iarg+1],"y") == 0) {
	cond[ncondition].lhs = Y;
	cond[ncondition].type = DOUBLE;
	cond[ncondition].stride = 3;
      } else if (strcmp(arg[iarg+1],"z") == 0) {
	cond[ncondition].lhs = Z;
	cond[ncondition].type = DOUBLE;
	cond[ncondition].stride = 3;
      } else if (arg[iarg+1][0] == 'i') {
	int index = atoi(&arg[iarg+1][1]);
	if (index < 1 || index > app->ninteger)
	  error->all(FLERR,
                     "Set if test on quantity application does not support");
	index--;
	cond[ncondition].lhs = IARRAY;
	cond[ncondition].type = INT;
	cond[ncondition].index = index;
	cond[ncondition].stride = 1;
      } else if (arg[iarg+1][0] == 'd') {
	int index = atoi(&arg[iarg+1][1]);
	if (index < 1 || index > app->ndouble)
	  error->all(FLERR,
                     "Set if test on quantity application does not support");
	index--;
	cond[ncondition].lhs = DARRAY;
	cond[ncondition].type = DOUBLE;
	cond[ncondition].index = index;
	cond[ncondition].stride = 1;
      } else error->all(FLERR,"Illegal set command");

      if (strcmp(arg[iarg+2],"<") == 0) cond[ncondition].op = LT;
      else if (strcmp(arg[iarg+2],"<=") == 0) cond[ncondition].op = LE;
      else if (strcmp(arg[iarg+2],">") == 0) cond[ncondition].op = GT;
      else if (strcmp(arg[iarg+2],">=") == 0) cond[ncondition].op = GE;
      else if (strcmp(arg[iarg+2],"=") == 0) cond[ncondition].op = EQ;
      else if (strcmp(arg[iarg+2],"!=") == 0) cond[ncondition].op = NEQ;
      else error->all(FLERR,"Illegal set command");

      if (cond[ncondition].type == INT)
	cond[ncondition].irhs = atoi(arg[iarg+3]);
      else cond[ncondition].drhs = atof(arg[iarg+3]);

      ncondition++;
      iarg += 4;

    } else error->all(FLERR,"Illegal set command");
  }

  if (domain->me == 0 && screen) fprintf(screen,"Setting site values ...\n");

  if (rhs == VALUE) set_single(lhs,rhs);
  else if (rhs == RANGE) set_range(lhs,rhs);
  else if (rhs == UNIQUE) set_single(lhs,rhs);
  else if (rhs == DISPLACE) set_displace(lhs,rhs);

  // statistics

  tagint nbig = count;
  tagint allcount;
  MPI_Allreduce(&nbig,&allcount,1,MPI_SPK_TAGINT,MPI_SUM,world);
    
  if (domain->me == 0) {
    if (screen) fprintf(screen,"  " TAGINT_FORMAT " settings made for %s\n",
			allcount,arg[0]);
    if (logfile) fprintf(logfile,"  " TAGINT_FORMAT " settings made for %s\n",
			 allcount,arg[0]);
  }

  memory->sfree(cond);
}

/* ----------------------------------------------------------------------
   set sites to a single value
   ivalue for VALUE, site ID for UNIQUE
   account for loop, region, fraction, condition options
------------------------------------------------------------------------- */

void Set::set_single(int lhs, int rhs)
{
  int i;
  tagint iglobal;

  int nlocal = app->nlocal;
  tagint maxID = app->max_site_ID();

  tagint *id = app->id;
  double **xyz = app->xyz;
  int **iarray = app->iarray;
  double **darray = app->darray;

  // if loopflag == 1, same RNG on every proc
  // if loopflag == 0, different RNG on every proc

  RandomPark *random = new RandomPark(ranmaster->uniform());
  if (loopflag == 0) {
    double seed = ranmaster->uniform();
    random->reset(seed,domain->me,100);
  }

  count = 0;

  if (loopflag) {
    std::map<tagint,int> hash;
    std::map<tagint,int>::iterator loc;
    for (i = 0; i < nlocal; i++) hash.insert(std::pair<tagint,int> (id[i],i));

    if (lhs == IARRAY) {
      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	  else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	    else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	  else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	    else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
	}
      }

    } else if (lhs == DARRAY) {
      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) darray[siteindex][i] = dvalue;
	  else darray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) darray[siteindex][i] = dvalue;
	    else darray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) darray[siteindex][i] = dvalue;
	  else darray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) darray[siteindex][i] = dvalue;
	    else darray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
	}
      }
    }

  } else {
    if (lhs == IARRAY) {
      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	  else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	    else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (random->uniform() >= fraction) continue;
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	  else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (random->uniform() >= fraction) continue;
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) iarray[siteindex][i] = ivalue;
	    else iarray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      }

    } else if (lhs == DARRAY) {
      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) darray[siteindex][i] = ivalue;
	  else darray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) darray[siteindex][i] = ivalue;
	    else darray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (random->uniform() < fraction) continue;
	  if (ncondition && condition(i)) continue;
	  if (rhs == VALUE) darray[siteindex][i] = ivalue;
	  else darray[siteindex][i] = id[i] % MAXSMALLINT;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (random->uniform() < fraction) continue;
	    if (ncondition && condition(i)) continue;
	    if (rhs == VALUE) darray[siteindex][i] = ivalue;
	    else darray[siteindex][i] = id[i] % MAXSMALLINT;
	    count++;
	  }
      }
    }
  }

  delete random;
}

/* ----------------------------------------------------------------------
   set sites to a range of values
   account for loop, region, fraction, condition options
------------------------------------------------------------------------- */

void Set::set_range(int lhs, int rhs)
{
  int i;
  tagint iglobal;

  int nlocal = app->nlocal;
  tagint maxID = app->max_site_ID();

  tagint *id = app->id;
  double **xyz = app->xyz;
  int **iarray = app->iarray;
  double **darray = app->darray;

  // if loopflag == 1, same RNG on every proc
  // if loopflag == 0, different RNG on every proc

  RandomPark *random = new RandomPark(ranmaster->uniform());
  if (loopflag == 0) {
    double seed = ranmaster->uniform();
    random->reset(seed,domain->me,100);
  }

  count = 0;

  if (loopflag) {
    std::map<tagint,int> hash;
    std::map<tagint,int>::iterator loc;
    for (i = 0; i < nlocal; i++) hash.insert(std::pair<tagint,int> (id[i],i));

    if (lhs == IARRAY) {
      int range = ivaluehi - ivaluelo + 1;

      if (regionflag == 0 && fraction == 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  ivalue = random->irandom(range);
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  iarray[siteindex][i] = ivalue-1 + ivaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  ivalue = random->irandom(range);
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    iarray[siteindex][i] = ivalue-1 + ivaluelo;
	    count++;
	  }
	}
      } else if (regionflag == 0 && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  ivalue = random->irandom(range);
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  iarray[siteindex][i] = ivalue-1 + ivaluelo;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  ivalue = random->irandom(range);
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    iarray[siteindex][i] = ivalue-1 + ivaluelo;
	    count++;
	  }
	}
      }

    } else if (lhs == DARRAY) {
      double range = dvaluehi - dvaluelo;

      if (regionflag == 0 && fraction == 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  dvalue = random->uniform();
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  darray[siteindex][i] = range*dvalue + dvaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  dvalue = random->uniform();
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    darray[siteindex][i] = range*dvalue + dvaluelo;
	    count++;
	  }
	}
      } else if (regionflag == 0 && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  dvalue = random->uniform();
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  darray[siteindex][i] = range*dvalue + dvaluelo;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (iglobal = 1; iglobal <= maxID; iglobal++) {
	  if (random->uniform() >= fraction) continue;
	  dvalue = random->uniform();
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    darray[siteindex][i] = range*dvalue + dvaluelo;
	    count++;
	  }
	}
      }
    }

  } else {
    if (lhs == IARRAY) {
      int range = ivaluehi - ivaluelo + 1;

      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  iarray[siteindex][i] = random->irandom(range)-1 + ivaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    iarray[siteindex][i] = random->irandom(range)-1 + ivaluelo;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (random->uniform() >= fraction) continue;
	  if (ncondition && condition(i)) continue;
	  iarray[siteindex][i] = random->irandom(range)-1 + ivaluelo;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (random->uniform() >= fraction) continue;
	    if (ncondition && condition(i)) continue;
	    iarray[siteindex][i] = random->irandom(range)-1 + ivaluelo;
	    count++;
	  }
      }

    } else if (lhs == DARRAY) {
      double range = dvaluehi - dvaluelo;

      if (regionflag == 0 && fraction == 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (ncondition && condition(i)) continue;
	  darray[siteindex][i] = range*random->uniform() + dvaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (ncondition && condition(i)) continue;
	    darray[siteindex][i] = range*random->uniform() + dvaluelo;
	    count++;
	  }
      } else if (regionflag == 0 && fraction < 1.0) {
	for (i = 0; i < nlocal; i++) {
	  if (random->uniform() < fraction) continue;
	  if (ncondition && condition(i)) continue;
	  darray[siteindex][i] = range*random->uniform() + dvaluelo;
	  count++;
	}
      } else if (regionflag && fraction < 1.0) {
	for (i = 0; i < nlocal; i++)
	  if (domain->regions[iregion]->match(xyz[i][0],xyz[i][1],xyz[i][2])) {
	    if (random->uniform() < fraction) continue;
	    if (ncondition && condition(i)) continue;
	    darray[siteindex][i] = range*random->uniform() + dvaluelo;
	    count++;
	  }
      }
    }
  }

  delete random;
}

/* ----------------------------------------------------------------------
   displace site coordinates
------------------------------------------------------------------------- */

void Set::set_displace(int lhs, int rhs)
{
}

/* ----------------------------------------------------------------------
   test local site I against if-test conditions
   return 0 if it meets all tests
   return 1 if it fails any if tests
------------------------------------------------------------------------- */

int Set::condition(int i)
{
  int *intarray;
  tagint *tagintarray;
  double *doublearray;

  for (int m = 0; m < ncondition; m++) {
    if (cond[m].lhs == IARRAY) intarray = app->iarray[cond[m].index];
    else if (cond[m].lhs == DARRAY) doublearray = app->darray[cond[m].index];
    else if (cond[m].lhs == X) doublearray = &app->xyz[0][0];
    else if (cond[m].lhs == Y) doublearray = &app->xyz[0][1];
    else if (cond[m].lhs == Z) doublearray = &app->xyz[0][2];
    else if (cond[m].lhs == ID) tagintarray = app->id;

    int offset = cond[m].stride * i;

    if (cond[m].op == LT) {
      if (cond[m].type == INT) {
	if (intarray[offset] >= cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] >= cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] >= cond[m].drhs) return 1;
      }
    } else if (cond[m].op == LE) {
      if (cond[m].type == INT) {
	if (intarray[offset] > cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] > cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] > cond[m].drhs) return 1;
      }
    } else if (cond[m].op == GT) {
      if (cond[m].type == INT) {
	if (intarray[offset] <= cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] <= cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] <= cond[m].drhs) return 1;
      }
    } else if (cond[m].op == GE) {
      if (cond[m].type == INT) {
	if (intarray[offset] < cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] < cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] < cond[m].drhs) return 1;
      }
    } else if (cond[m].op == EQ) {
      if (cond[m].type == INT) {
	if (intarray[offset] != cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] != cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] != cond[m].drhs) return 1;
      }
    } else if (cond[m].op == NEQ) {
      if (cond[m].type == INT) {
	if (intarray[offset] == cond[m].irhs) return 1;
      } else if (cond[m].type == TAGINT) {
	if (tagintarray[offset] == cond[m].irhs) return 1;
      } else {
	if (doublearray[offset] == cond[m].drhs) return 1;
      }
    }
  }
  return 0;
}
