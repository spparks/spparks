/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
 ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_surf.h"
#include "solve.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS;

enum{DEPOSITION,HOPLEFT,HOPRIGHT};

#define DELTA 10
#define TOL 1.0e-5

#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

AppSurf::AppSurf(SPK *spk, int narg, char **arg) : App(spk,narg,arg)
{
  if (narg != 4) error->all("Invalid app_style surf command");

  nlattice = atoi(arg[1]);
  strain = atof(arg[2]);
  seed = atoi(arg[3]);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nlocal = nghost = maxatom = 0;
  atoms = NULL;

  nevents = maxevent = 0;
  events = NULL;
  rates = NULL;

  random = new RandomPark(seed);
  attempt_frequency = 1.0;

  epsilon = sigma = 1.0;
  sigma6 = sigma12 = 1.0;
   
  // default settings

  ntimestep = 0;
  nstats = ndump = 0;
  temperature = 0.0;
  cutoff = 2.5;
  rate_deposit = 1.0;
  fp = NULL;
}

/* ---------------------------------------------------------------------- */

AppSurf::~AppSurf()
{
  memory->sfree(atoms);
  memory->sfree(events);
  memory->sfree(rates);

  delete random;
  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void AppSurf::init()
{
  // setup initial lattice
  // # of substrate layers depends on cutoff

  double hspacing = strain * pow(2.0,1.0/6.0);
  double vspacing = sqrt(3.0)/2.0 * hspacing;
  dist_hop = hspacing;

  cutsq = cutoff * cutoff;
  double cutneigh = 1.2 * hspacing;
  cutneighsq = cutneigh * cutneigh;

  xlo = 0.0;
  xhi = nlattice * hspacing;
  xprd = xhi - xlo;

  int nlayers = static_cast<int> (cutoff/vspacing) + 1;

  int i,j;
  double xtmp,ztmp;

  for (j = 0; j < nlayers; j++) {
    for (i = 0; i < nlattice; i++) {
      xtmp = i*hspacing + (j % 2)*(0.5*hspacing);
      ztmp = j*vspacing;
      add_atom(nlocal,nlocal+1,1,xtmp,ztmp);
      nlocal++;
    }
  }

  zlo = 0.0;
  zhi = zmax();

  // setup future stat and dump calls

  stats_next = nstats;
  dump_next = ndump;

  // print stats header
  
  if (me == 0) {
    if (screen) {
      fprintf(screen,"Timestep Time Natoms Energy");
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"Timestep Time Natoms Energy");
      fprintf(logfile,"\n");
    }
  }

  // initial stats and dump

  energy = relax();
  if (nstats) stats();
  if (ndump) dump();
}

/* ---------------------------------------------------------------------- */

void AppSurf::input(char *command, int narg, char **arg)
{
  if (narg == 0) error->all("Invalid command");
  if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"potential") == 0) set_potential(narg,arg);
  else if (strcmp(command,"rates") == 0) set_rates(narg,arg);
  else if (strcmp(command,"stats") == 0) set_stats(narg,arg);
  else if (strcmp(command,"dump") == 0) set_dump(narg,arg);
  else error->all("Invalid command");
}

/* ----------------------------------------------------------------------
   perform a run
 ------------------------------------------------------------------------- */

void AppSurf::run(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  nsteps = atoi(arg[0]);
  time = 0.0;

  // init classes used by this app
  
  init();
  timer->init();
  
  // perform the run
  
  iterate();
  
  // final statistics
  
  Finish finish(spk);
}

/* ----------------------------------------------------------------------
   iterate on events
 ------------------------------------------------------------------------- */

void AppSurf::iterate()
{
  int i,which;
  double delta,rate,dt,xtmp,eng,fx,fz;

  timer->barrier_start(TIME_LOOP);

  for (int m = 0; m < nsteps; m++) {
    ntimestep++;

    ghost_comm();

    nevents = 0;
    add_event(nlocal,DEPOSITION,rate_deposit);

    for (i = 0; i < nlocal; i++) {
      if (atoms[i].type == 1) continue;
      else if (atoms[i].type == 3) atoms[i].type == 2;

      if (count_neigh(i) >= 4) continue;
      
      delta = find_barrier(i,atoms[i].x - dist_hop);
      rate = attempt_frequency * exp(-delta / temperature);
      add_event(i,HOPLEFT,rate);

      delta = find_barrier(i,atoms[i].x + dist_hop);
      rate = attempt_frequency * exp(-delta / temperature);
      add_event(i,HOPRIGHT,rate);
    }

    solve->init(nevents,rates);
    which = solve->event(&dt);
    time += dt;

    if (events[which].style == DEPOSITION) {
      xtmp = xlo + random->uniform()*xprd;
      add_atom(nlocal,nlocal+1,3,xtmp,0.0);
      nlocal++;

      ghost_comm();
      atoms[nlocal-1].z = zrelax(nlocal-1,eng,fx,fz);

      double ztmp = atoms[nlocal-1].z;
      atoms[nlocal-1].z = ztmp - sigma/100.0;
      double fzlo;
      double elo = engforce(nlocal-1,fx,fzlo);
      atoms[nlocal-1].z = ztmp + sigma/100.0;
      double fzhi;
      double ehi = engforce(nlocal-1,fx,fzhi);
      printf("DEP: lo = %g,%g; mid = %g,%g; hi = %g,%g\n",
	     elo,fzlo,eng,fz,ehi,fzhi);
      atoms[nlocal-1].z = ztmp;

    } else {
      i = events[which].iatom;
      if (events[which].style == HOPLEFT) atoms[i].x -= dist_hop;
      else atoms[i].x += dist_hop;

      if (atoms[i].x < xlo) atoms[i].x += xprd;
      else if (atoms[i].x >= xhi) atoms[i].x -= xprd;
      ghost_comm();
      atoms[i].z = zrelax(i,eng,fx,fz);
    }

    energy = relax();
    pbc();

    if (ntimestep == stats_next) {
      timer->stamp();
      stats();
      stats_next += nstats;
      timer->stamp(TIME_OUTPUT);
    }
    if (ntimestep == dump_next) {
      timer->stamp();
      dump();
      dump_next += ndump;
      timer->stamp(TIME_OUTPUT);
    }
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   count # of neighs for IATOM within neighbor cutoff
------------------------------------------------------------------------- */

int AppSurf::count_neigh(int iatom)
{
  double dx,dz,rsq;

  int count = 0;
  for (int i = 0; i < nlocal+nghost; i++) {
    if (i == iatom) continue;
    dx = atoms[i].x - atoms[iatom].x;
    dz = atoms[i].z - atoms[iatom].z;
    rsq = dx*dx + dz*dz;
    if (rsq < cutneighsq) count++;
  }
  return count;
}

/* ----------------------------------------------------------------------
   find highest z of all owned atoms
------------------------------------------------------------------------- */

double AppSurf::zmax()
{
  zhi = 0.0;
  for (int i = 0; i < nlocal; i++)
    zhi = MAX(zhi,atoms[i].z);
  return zhi;
}

/* ----------------------------------------------------------------------
   relax Z position of IATOM while holding all other atoms fixed
   return new Z coord
   also compute final eng, fx, fz of IATOM
------------------------------------------------------------------------- */

double AppSurf::zrelax(int iatom, double &eng, double &fx, double &fz)
{
  double z_initial = atoms[iatom].z;

  // dbrent requires bracketed input:
  // lo < mid < hi, eng(mid) < eng(lo), eng(mid) < eng(hi)
  // strategy:
  //   start z at zmax + sigma
  //   reduce z by delta until gradient is up, this is lo
  //   increase z by delta until gradient is down, this is hi
  //   whichever elo/ehi is smaller becomes mid, epsilon away is new lo/hi

  double z,lo,mid,hi,elo,emid,ehi;

  double delta = sigma/2.0;
  double eps = sigma/100.0;

  z = zmax() + delta;
  atoms[iatom].z = z;
  elo = engforce(iatom,fx,fz);
  while (fz <= 0.0) {
    z -= delta;
    atoms[iatom].z = z;
    elo = engforce(iatom,fx,fz);
  }
  lo = z;

  z += delta;
  atoms[iatom].z = z;
  ehi = engforce(iatom,fx,fz);
  while (fz >= 0.0) {
    z += delta;
    atoms[iatom].z = z;
    ehi = engforce(iatom,fx,fz);
  }
  hi = z;

  if (elo < ehi) {
    mid = lo;
    lo = mid - eps;
  } else {
    mid = hi;
    hi = mid + eps;
  }

  // dbrent will compute optimal znew
  // compute eng,fx,fz at znew to return to caller

  double znew;
  double tmp = dbrent(iatom,lo,mid,hi,TOL,&znew);
  atoms[iatom].z = znew;
  eng = engforce(iatom,fx,fz);
  
  atoms[iatom].z = z_initial;
  return znew;
}

/* ----------------------------------------------------------------------
   add an atom to atom list at location IATOM
------------------------------------------------------------------------- */

void AppSurf::add_atom(int iatom, int id, int type, double x, double z)
{
  if (iatom == maxatom) {
    maxatom += DELTA;
    atoms = (OneAtom *) 
      memory->srealloc(atoms,maxatom*sizeof(OneAtom),"surf:atoms");
  }

  atoms[iatom].id = id;
  atoms[iatom].type = type;
  atoms[iatom].x = x;
  atoms[iatom].z = z;
}

/* ----------------------------------------------------------------------
   add an event to event list
------------------------------------------------------------------------- */

void AppSurf::add_event(int iatom, int style, double rate)
{
  if (nevents == maxevent) {
    maxevent += DELTA;
    events = (OneEvent *) 
      memory->srealloc(events,maxevent*sizeof(OneEvent),"surf:events");
    rates = (double *) 
      memory->srealloc(rates,maxevent*sizeof(double),"surf:rates");
  }

  events[nevents].iatom = iatom;
  events[nevents].style = style;
  rates[nevents] = rate;
  nevents++;
}

/* ----------------------------------------------------------------------
   find a diffusive barrier for IATOM hopping to xnew
   return barrier height
------------------------------------------------------------------------- */

double AppSurf::find_barrier(int iatom, double xnew)
{
  double eng_initial,x_initial,z_initial;
  double lo,hi,mid;
  double znew,eng,fx,fz;

  eng_initial = relax();

  x_initial = atoms[iatom].x;
  z_initial = atoms[iatom].z;

  lo = MIN(x_initial,xnew);
  hi = MAX(x_initial,xnew);

  int niterate = 0;
  while (niterate < 10) {
    mid = 0.5 * (lo+hi);
    atoms[iatom].x = mid;
    znew = zrelax(iatom,eng,fx,fz);
    
    if (fx < 0.0) lo = mid;
    else if (fx > 0.0) hi = mid;
    else break;

    niterate++;
  }

  atoms[iatom].x = x_initial;
  atoms[iatom].z = z_initial;

  return eng - eng_initial;
}

/* ----------------------------------------------------------------------
   relax the entire system
   return total energy of system
------------------------------------------------------------------------- */

double AppSurf::relax()
{
  int i,j;
  double dx,dz,rsq,r2inv,r6inv;

  double lj3 = 4.0 * epsilon * pow(sigma,12.0);
  double lj4 = 4.0 * epsilon * pow(sigma,6.0);
  double ratio = sigma / cutoff;
  double offset = 4.0 * epsilon * (pow(ratio,12.0) - pow(ratio,6.0));

  double eng = 0.0;
  for (i = 0; i < nlocal; i++) {
    for (j = i+1; j < nlocal+nghost; j++) {
      if (atoms[i].type == 1 && atoms[j].type == 1) continue;
      dx = atoms[i].x - atoms[j].x;
      dz = atoms[i].z - atoms[j].z;
      rsq = dx*dx + dz*dz;
      if (rsq >= cutsq) continue;
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;
      eng += r6inv * (lj3*r6inv - lj4) - offset;
    }
  }

  return eng;
}

/* ----------------------------------------------------------------------
   compute eng and force on IATOM
   return eng, set fx & fz
------------------------------------------------------------------------- */

double AppSurf::engforce(int iatom, double &fx, double &fz)
{
  double dx,dz,rsq,r2inv,r6inv,fforce;

  double lj1 = 48.0 * epsilon * pow(sigma,12.0);
  double lj2 = 24.0 * epsilon * pow(sigma,6.0);
  double lj3 = 4.0 * epsilon * pow(sigma,12.0);
  double lj4 = 4.0 * epsilon * pow(sigma,6.0);
  double ratio = sigma / cutoff;
  double offset = 4.0 * epsilon * (pow(ratio,12.0) - pow(ratio,6.0));

  fx = fz = 0.0;
  double eng = 0.0;

  for (int j = 0; j < nlocal+nghost; j++) {
    if (j == iatom) continue;
    dx = atoms[iatom].x - atoms[j].x;
    dz = atoms[iatom].z - atoms[j].z;
    rsq = dx*dx + dz*dz;
    if (rsq >= cutsq) continue;
    r2inv = 1.0/rsq;
    r6inv = r2inv*r2inv*r2inv;
    fforce = r6inv * (lj1*r6inv - lj2) * r2inv;
    fx += dx * fforce;
    fz += dz * fforce;
    eng += r6inv * (lj3*r6inv - lj4) - offset;
  }

  return eng;
}

/* ----------------------------------------------------------------------
   acquire ghost atoms (beyond PBC on one proc)
------------------------------------------------------------------------- */

void AppSurf::ghost_comm()
{
  nghost = 0;
  for (int i = 0; i < nlocal; i++) {
    if (atoms[i].x < xlo+cutoff) {
      add_atom(nlocal+nghost,atoms[i].id,atoms[i].type,
	       atoms[i].x+xprd,atoms[i].z);
      nghost++;
    } else if (atoms[i].x > xhi-cutoff) {
      add_atom(nlocal+nghost,atoms[i].id,atoms[i].type,
	       atoms[i].x-xprd,atoms[i].z);
      nghost++;
    }
  }
}

/* ----------------------------------------------------------------------
   enforce PBC for owned atoms
------------------------------------------------------------------------- */

void AppSurf::pbc()
{
  for (int i = 0; i < nlocal; i++) {
    if (atoms[i].x < xlo) atoms[i].x += xprd;
    else if (atoms[i].x >= xhi) atoms[i].x -= xprd;
  }
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppSurf::stats()
{
  if (me == 0) {
    if (screen) {
      fprintf(screen,"%d %g %d %g\n",ntimestep,time,nlocal,energy);
    }
    if (logfile) {
      fprintf(logfile,"%d %g %d %g\n",ntimestep,time,nlocal,energy);
    }
  }
}

/* ----------------------------------------------------------------------
   dump a snapshot of atoms
------------------------------------------------------------------------- */

void AppSurf::dump()
{
  zhi = zmax();

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",nlocal);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",xlo,xhi);
  fprintf(fp,"%g %g\n",zlo,zhi);
  fprintf(fp,"%g %g\n",-0.5,0.5);
  fprintf(fp,"ITEM: ATOMS\n");

  // write atoms

  for (int i = 0; i < nlocal; i++)
    fprintf(fp,"%d %d %g %g %g\n",
	    atoms[i].id,atoms[i].type,atoms[i].x,atoms[i].z,0.0);
}

/* ---------------------------------------------------------------------- */

void AppSurf::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal temperature command");
  temperature = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppSurf::set_stats(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal stats command");
  nstats = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppSurf::set_potential(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal potential command");
  cutoff = atof(arg[1]);
}

/* ---------------------------------------------------------------------- */

void AppSurf::set_rates(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal rates command");
  rate_deposit = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void AppSurf::set_dump(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal dump command");
  ndump = atoi(arg[0]);
  if (me == 0) {
    fp = fopen(arg[1],"w");
    if (fp == NULL) error->one("Cannot open dump file");
  }
}

/* ----------------------------------------------------------------------
   1d minimization with gradient info
   slightly modified dbrent() routine from Numerical Recipes
   alters the z coord of IATOM
   return new min eng
   returns xmin = new z coord where eng is min
------------------------------------------------------------------------- */

double AppSurf::dbrent(int iatom, double ax, double bx, double cx,
		       double tol, double *xmin)
{
  int ok1,ok2;
  double a,b,d,d1,d2,du,dv,dw,dx,e;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  double xforce,zforce;

  e = 0.0;
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;

  atoms[iatom].z = x;
  fw = fv = fx = engforce(iatom,xforce,zforce);
  dw = dv = dx = zforce;

  for (int iter = 0; iter < ITMAX; iter++) {
    xm = 0.5*(a+b);
    tol1 = tol*fabs(x) + ZEPS;
    tol2 = 2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1 = 2.0*(b-a);
      d2 = d1;
      if (dw != dx) d1 = (w-x)*dx/(dx-dw);
      if (dv != dx) d2 = (v-x)*dx/(dx-dv);
      u1 = x+d1;
      u2 = x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde = e;
      e = d;
      if (ok1 || ok2) {
	if (ok1 && ok2)
	  d = (fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d = d1;
	else
	  d = d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u = x+d;
	  if (u-a < tol2 || b-u < tol2) d = SIGN(tol1,xm-x);
	} else {
	  d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
    }

    if (fabs(d) >= tol1) {
      u = x + d;
      atoms[iatom].z = u;
      fu = engforce(iatom,xforce,zforce);
    } else {
      u = x + SIGN(tol1,d);
      atoms[iatom].z = u;
      fu = engforce(iatom,xforce,zforce);
      if (fu > fx) {
	*xmin = x;
	return fx;
      }
    }
    du = zforce;

    if (fu <= fx) {
      if (u >= x) a = x;
      else b = x;
      MOV3(v,fv,dv, w,fw,dw)
      MOV3(w,fw,dw, x,fx,dx)
      MOV3(x,fx,dx, u,fu,du)
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw)
	MOV3(w,fw,dw, u,fu,du)
      } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du)
      }
    }
  }

  error->warning("Too many iterations in routine dbrent");
  return 0.0;
}
