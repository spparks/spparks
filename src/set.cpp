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
#include "stdio.h"
#include "unistd.h"
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

#ifdef SPK_STITCH
#include "stitch.h"
#endif

using namespace SPPARKS_NS;

enum{SITE,IARRAY,DARRAY,X,Y,Z,XYZ,ID};
enum{VALUE,RANGE,UNIQUE,DISPLACE,STITCH,BFILE};
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
  filename = NULL;
  tstamp = NULL;

  if (strcmp(arg[0],"site") == 0) {
    lhs = SITE;
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
    if (lhs == SITE || lhs == IARRAY) ivalue = atoi(arg[2]);
    else if (lhs == DARRAY) dvalue = atof(arg[2]);
    else error->all(FLERR,"Illegal set command");
    iarg = 3;
  } else if (strcmp(arg[1],"range") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal set command");
    rhs = RANGE;
    if (lhs == SITE || lhs == IARRAY) {
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
    
#ifdef SPK_STITCH
  } else if (strcmp(arg[1],"stitch") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal set command");
    rhs = STITCH;
    int n = strlen(arg[2]) + 1;
    filename = new char[n];
    strcpy(filename,arg[2]);
    n = strlen(arg[3]) + 1;
    tstamp = new char[n];
    strcpy(tstamp,arg[3]);
    iarg = 4;
#endif
    
  } else if (strcmp(arg[1],"bfile") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal set command");
    rhs = BFILE;
    int n = strlen(arg[2]) + 1;
    filename = new char[n];
    strcpy(filename,arg[2]);
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
      cond = (Condition *)
	memory->srealloc(cond, (ncondition+1)*sizeof(Condition), "set:cond");
      if (strcmp(arg[iarg+1],"id") == 0) {
         cond[ncondition].lhs = ID;
         cond[ncondition].type = TAGINT;
         cond[ncondition].stride = 1;
      } else if (strcmp(arg[iarg+1],"site") == 0) {
         cond[ncondition].lhs = IARRAY;
         cond[ncondition].type = INT;
         cond[ncondition].index = 0;
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

      if (cond[ncondition].type == INT) cond[ncondition].irhs = atoi(arg[iarg+3]);
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
  else if (rhs == STITCH) set_stitch(lhs,rhs);
  else if (rhs == BFILE) set_binary_file(lhs,rhs);

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

  // cleanup

  delete [] filename;
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
  tagint minID = app->min_site_ID();
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
    MyHash hash;
    MyIterator loc;

    for (i = 0; i < nlocal; i++) hash.insert(std::pair<tagint,int> (id[i],i));

    if (lhs == SITE || lhs == IARRAY) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
    if (lhs == SITE || lhs == IARRAY) {
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
  tagint minID = app->min_site_ID();
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
    MyHash hash;
    MyIterator loc;

    for (i = 0; i < nlocal; i++) hash.insert(std::pair<tagint,int> (id[i],i));

    if (lhs == SITE || lhs == IARRAY) {
      int range = ivaluehi - ivaluelo + 1;

      if (regionflag == 0 && fraction == 1.0) {
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
	  ivalue = random->irandom(range);
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  iarray[siteindex][i] = ivalue-1 + ivaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
	  dvalue = random->uniform();
	  loc = hash.find(iglobal);
	  if (loc == hash.end()) continue;
	  i = loc->second;
	  if (ncondition && condition(i)) continue;
	  darray[siteindex][i] = range*dvalue + dvaluelo;
	  count++;
	}
      } else if (regionflag && fraction == 1.0) {
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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
	for (iglobal = minID; iglobal <= maxID; iglobal++) {
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

    if (lhs == SITE || lhs == IARRAY) {
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
   displace site coordinates - NOT YET implemented
------------------------------------------------------------------------- */

void Set::set_displace(int lhs, int rhs) {}

/* ----------------------------------------------------------------------
   set sites from values in stitch file
------------------------------------------------------------------------- */

void Set::set_stitch(int lhs, int rhs)
{
#ifdef SPK_STITCH

  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Set stitch only allowed for on-lattice apps");

  applattice = (AppLattice *) app;

  if (!applattice->simple)
    error->all(FLERR,"Set stitch requires simple square or cubic lattice");
  if (lhs != SITE && lhs != IARRAY)
    error->all(FLERR,"Set stitch only supports integer values");

  StitchFile *stitch_file;
  int err = stitch_open(&stitch_file,MPI_COMM_WORLD,filename);
  // TODO: process err
  //

  double time;
  if (strcmp(tstamp,"first") == 0 || strcmp(tstamp,"last") == 0) {
    double first_time,last_time,tmp;
    int64_t int_tmp=-1;
    err = stitch_get_parameters(stitch_file,&tmp,&tmp,&int_tmp,&first_time,&last_time);
     // TODO: process err
     //
    if (strcmp(tstamp,"first") == 0) time = first_time;
    else time = last_time;
  } else time = atof(tstamp);

  int block[6], uninitialized_bb[6];
  block[0] = applattice->xlo_me_simple;
  block[1] = applattice->xhi_me_simple+1;
  block[2] = applattice->ylo_me_simple;
  block[3] = applattice->yhi_me_simple+1;
  block[4] = applattice->zlo_me_simple;
  block[5] = applattice->zhi_me_simple+1;
  int xlo = applattice->xlo_me_simple;
  int xhi = applattice->xhi_me_simple+1;
  int ylo = applattice->ylo_me_simple;
  int yhi = applattice->yhi_me_simple+1;
  int zlo = applattice->zlo_me_simple;
  int zhi = applattice->zhi_me_simple+1;

  // get field id
  
  int64_t field_id;
  bool init=false;
  {
     int err;
     char label[8];
     if (SITE==lhs) {
        sprintf(label,"site");
     } else if (IARRAY==lhs) {
        sprintf(label,"i%d",siteindex+1);
     } else if (DARRAY==lhs) {
        sprintf(label,"d%d",siteindex+1);
     } else {
       error->all(FLERR,"Set stitch command failed due to "
		  "invalid field specification; "
		  "must be 'site', i1,i2,..., or d1,d2,d3,....");
     }
     err=stitch_query_field (stitch_file,label,&field_id);
     if (STITCH_FIELD_NOT_FOUND_ID==field_id) {
        // then create field 
        if (SITE==lhs || IARRAY==lhs) {
           // Specify value type of value and its undefined value
           union StitchTypesUnion v;v.i32=-1;
           // Scalar value
           int32_t scalar=1;
           err = stitch_create_field(stitch_file,label,STITCH_INT32,
                 v,scalar,&field_id);
           // TODO: process err
           //
        } else if (DARRAY==lhs) {
           union StitchTypesUnion v;v.f64=-1;
           // Scalar value
           int32_t scalar=1;
           err = stitch_create_field (stitch_file,label,STITCH_FLOAT64,
			    v,scalar,&field_id);
           // TODO: process err
           //
        } else {
          error->all(FLERR,"'set stitch' command failed on read; "
           "Specified field does not exist; also unable to create field type specified.");
        }
        // field was initialized on 'stitch' file; set 'init'=true
        init=true;
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        if (0==my_rank){
           printf("'set stitch' creating field '%s' on 'stitch file'\n",label);
        }
     }
  }

  {
     int err;

     // get field data
     int32_t is_new_time=-1;
     if (lhs == SITE) {
        int32_t *idata = app->iarray[0];
        err =stitch_read_block_int32(stitch_file,field_id,&time,
				     block,idata,&is_new_time);
     } else if (lhs == IARRAY) {
        int32_t *idata = app->iarray[siteindex];
        //For read:
        // read_flag == True when time does not exist in the file,
	//   but data may be returned based on older times
        // read_flag == False when exists in the database. The actual data returned
	//   will vary based on what is selected. It could be nothing.
        // int stitch_read_block_int32 (const StitchFile * file, int64_t field_id,
	//   double * time, int32_t * bb, int32_t * buffer, int32_t * is_new_time);
        err = stitch_read_block_int32(stitch_file,field_id,&time,
				      block,idata,&is_new_time);
     } else if (lhs == DARRAY) {
        double *real_data = app->darray[siteindex];
        err = stitch_read_block_float64(stitch_file,field_id,&time,
					block,real_data,&is_new_time);
     }
     // TODO: process err
     // if 'field' was newly created above, then 'init' is true 
     //    then site values 'spin' values will have default values 
     //    set by stitch library; default values must be handled by 
     //    app
     if (is_new_time && false==init)
       error->all(FLERR,"Set stitch command failed on read_block; "
		  "libstitch did not match a time.");

     err = stitch_close(&stitch_file);
     // TODO: process err
     //
  }

  /*
  int read_flag=-1;
  printf("0:Stitch read_block; time=%3.1f, read_flag=%d\n",time,read_flag);
  err = stitch_read_block(stitch_file,&time,block,idata,uninitialized_bb,&read_flag);
  printf("1:Stitch read_block; time=%3.1f, read_flag=%d\n",time,read_flag);
  printf("2: block xlo, xhi, ylo, yhi, zlo, zli = %5d,%5d,%5d,%5d,%5d,%5d",xlo,xhi,ylo,yhi,zlo,zhi);
  if(read_flag)
    error->all(FLERR,"Set stitch command failed on read; libstitch did not match a time.");
  err = stitch_close(&stitch_file);
  */
  
  count = app->nlocal;

#endif
}

/* ----------------------------------------------------------------------
   set sites from values in a binary file
------------------------------------------------------------------------- */

void Set::set_binary_file(int lhs, int rhs)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Set bfile only allowed for on-lattice apps");
 
  applattice = (AppLattice *) app;

  if (!applattice->simple)
    error->all(FLERR,"Set bfile requires simple square or cubic lattice");

  int me;
  MPI_Comm_rank(world,&me);

  // xlo to zhi = extent of lattice indices for entire box
  // xlo_me to zhi_me = extent of lattice indices for my subdomain
  // nxyz = size of simulation box lattice
  // nxyz_me = size of my portion of simulation box lattice

  int xlo,xhi,ylo,yhi,zlo,zhi;
  xlo = applattice->xlo_simple;
  xhi = applattice->xhi_simple;
  ylo = applattice->ylo_simple;
  yhi = applattice->yhi_simple;
  zlo = applattice->zlo_simple;
  zhi = applattice->zhi_simple;

  int xlo_me,xhi_me,ylo_me,yhi_me,zlo_me,zhi_me;
  xlo_me = applattice->xlo_me_simple;
  xhi_me = applattice->xhi_me_simple;
  ylo_me = applattice->ylo_me_simple;
  yhi_me = applattice->yhi_me_simple;
  zlo_me = applattice->zlo_me_simple;
  zhi_me = applattice->zhi_me_simple;

  int nx = xhi-xlo + 1;
  int ny = yhi-ylo + 1;
  int nz = zhi-zlo + 1;
  int nx_me = xhi_me-xlo_me + 1;
  int ny_me = yhi_me-ylo_me + 1;
  int nz_me = zhi_me-zlo_me + 1;

  // read header of binary file
  // insure it matches current simulation box

  int i,m,ix,iy,iz,xfile,yfile,zfile;
  int nxyz[3];
  FILE *fp;

  if (me == 0) {
    fp = fopen(filename,"rb");
    if (fp == NULL) error->one(FLERR,"Set could not open binary file");
    fread(nxyz,sizeof(int),3,fp);
    if (nxyz[0] != nx || nxyz[1] != ny || nxyz[2] != nz)
      error->one(FLERR,"Set binary file does not match current lattice");
  }

  // proc 0 reads all data from file, bcasts it to all procs
  // each proc extracts the subset for its sub-domain sites
  // assume file stores info for all sites: x fastest, then y, z slowest

  int nlocal = app->nlocal;

  if (lhs == SITE || lhs == IARRAY) {
    int *buf;
    memory->create(buf,nx*ny*nz,"set:buf");

    if (me == 0) fread(buf,sizeof(int),nx*ny*nz,fp);
    MPI_Bcast(buf,nx*ny*nz,MPI_INT,0,world);

    int *array = app->iarray[siteindex];
    for (i = 0; i < nlocal; i++) {
      ix = i % nx_me;
      iy = (i/nx_me) % ny_me;
      iz = i / (nx_me*ny_me);
      xfile = xlo_me-xlo + ix;
      yfile = ylo_me-ylo + iy;
      zfile = zlo_me-zlo + iz;
      m = nx*ny*zfile + nx*yfile + xfile;
      array[i] = buf[m];
    }

    memory->destroy(buf);

  } else if (lhs == DARRAY) {
    double *buf;
    memory->create(buf,nx*ny*nz,"set:buf");

    if (me == 0) fread(buf,sizeof(double),nx*ny*nz,fp);
    MPI_Bcast(buf,nx*ny*nz,MPI_INT,0,world);

    double *array = app->darray[siteindex];
    for (i = 0; i < nlocal; i++) {
      ix = i % nx_me;
      iy = (i/nx_me) % ny_me;
      iz = i / (nx_me*ny_me);
      xfile = xlo_me-xlo + ix;
      yfile = ylo_me-ylo + iy;
      zfile = zlo_me-zlo + iz;
      m = nx*ny*zfile + nx*yfile + xfile;
      array[i] = buf[m];
    }

    memory->destroy(buf);
  }

  if (me == 0) fclose(fp);

  count = nlocal;
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
