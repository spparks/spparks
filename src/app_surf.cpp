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

#define DELTA_ATOM  10
#define DELTA_EVENT 2
#define DELTA_NEIGH 100
#define TOL 1.0e-5

#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MAX_ITER 1000
#define MAX_EVAL 10000
#define MAX_LINE 10
#define DMIN 0.001
#define DMAX 0.1
#define TOLERANCE 1.0e-6

#define EPS         1.0e-10
#define SCAN_FACTOR 2.0
#define SECANT_EPS  1.0e-6

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//#define DEBUG 1

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

  maxatomneigh = maxneigh = 0;
  firstneigh = NULL;
  numneigh = NULL;
  neigh = NULL;

  nevents = maxevent = 0;
  events = NULL;
  rates = NULL;

  random = new RandomPark(seed);
  attempt_frequency = 1.0;

  epsilon = sigma = 1.0;

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
  delete [] firstneigh;
  delete [] numneigh;
  memory->sfree(neigh);
  memory->sfree(events);
  memory->sfree(rates);

  delete random;
  if (fp) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void AppSurf::init()
{
  lj1 = 48.0 * epsilon * pow(sigma,12.0);
  lj2 = 24.0 * epsilon * pow(sigma,6.0);
  lj3 = 4.0 * epsilon * pow(sigma,12.0);
  lj4 = 4.0 * epsilon * pow(sigma,6.0);
  double ratio = sigma / cutoff;
  offset = 4.0 * epsilon * (pow(ratio,12.0) - pow(ratio,6.0));

  // cutoffs

  cutsq = cutoff * cutoff;
  skin = sigma;
  cutneigh = cutoff + skin;
  cutneighsq = cutneigh * cutneigh;

  // setup initial lattice
  // # of substrate layers depends on cutoff

  double hspacing = strain * pow(2.0,1.0/6.0);
  double vspacing = sqrt(3.0)/2.0 * hspacing;
  hop_distance = hspacing;

  cutcount = 1.2 * hspacing;
  cutcountsq = cutcount * cutcount;

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

  zlo = zhi = 0.0;
  for (i = 0; i < nlocal; i++) zhi = MAX(zhi,atoms[i].z);

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

  energy = etotal();
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
  int i,iatom,which;
  double delta,rate,dt,xnew,znew;
  double eng,fx,fz;
  double eps = sigma/100.0;

  timer->barrier_start(TIME_LOOP);

  for (int m = 0; m < nsteps; m++) {
    ntimestep++;

    ghost_comm();

    nevents = 0;
    add_event(nlocal,DEPOSITION,rate_deposit,0.0,0.0);

    for (i = 0; i < nlocal; i++) {
      if (atoms[i].type == 1) continue;
      else if (atoms[i].type == 3) atoms[i].type = 2;

      if (count_neigh(i) >= 4) continue;

      neighbor_z(i,hop_distance);

      delta = find_barrier(i,atoms[i].x - hop_distance,&xnew,&znew);
      rate = attempt_frequency * exp(-delta / temperature);
      add_event(i,HOPLEFT,rate,xnew-eps,znew);

      delta = find_barrier(i,atoms[i].x + hop_distance,&xnew,&znew);
      rate = attempt_frequency * exp(-delta / temperature);
      add_event(i,HOPRIGHT,rate,xnew+eps,znew);
    }

    solve->init(nevents,rates);
    which = solve->event(&dt);
    time += dt;

    if (events[which].style == DEPOSITION) {
      iatom = nlocal;
      xnew = xlo + random->uniform()*xprd;
      add_atom(nlocal,nlocal+1,3,xnew,0.0);
      nlocal++;

      ghost_comm();
      neighbor_z(iatom,0.0);
      atoms[iatom].z = zrelax(iatom,&eng,&fx,&fz);

#ifdef DEBUG
      double ztmp = atoms[iatom].z;
      atoms[iatom].z = ztmp - sigma/100.0;
      double elo = eng_force(1,&iatom);
      double fzlo = atoms[iatom].fz;
      atoms[iatom].z = ztmp + sigma/100.0;
      double ehi = eng_force(1,&iatom);
      double fzhi = atoms[iatom].fz;
      printf("DEP: lo = %g,%g; mid = %g,%g; hi = %g,%g\n",
	     elo,fzlo,eng,fz,ehi,fzhi);
      atoms[iatom].z = ztmp;
#endif

    } else {
      iatom = events[which].iatom;
      atoms[iatom].x = events[which].x;
      atoms[iatom].z = events[which].z;

      if (atoms[iatom].x < xlo) atoms[iatom].x += xprd;
      else if (atoms[iatom].x >= xhi) atoms[iatom].x -= xprd;
      ghost_comm();

#ifdef DEBUG
      double xtmp = atoms[iatom].x;
      double ztmp = atoms[iatom].z;
      double xmid = xtmp;
      double zmid = ztmp;
      if (events[which].style == HOPLEFT) xmid += eps;
      if (events[which].style == HOPRIGHT) xmid -= eps;
      neighbor_z(iatom,eps);

      atoms[iatom].x = xmid;
      atoms[iatom].z = zmid;
      double emid = eng_force(1,&iatom);
      double fxmid = atoms[iatom].fx;
      double fzmid = atoms[iatom].fz;

      atoms[iatom].x = xmid;
      atoms[iatom].z = zmid - sigma/100.0;
      double ezlo = eng_force(1,&iatom);
      double fzlo = atoms[iatom].fz;

      atoms[iatom].x = xmid;
      atoms[iatom].z = zmid + sigma/100.0;
      double ezhi = eng_force(1,&iatom);
      double fzhi = atoms[iatom].fz;

      atoms[iatom].x = xmid - sigma/100.0;
      atoms[iatom].z = zmid;
      double exlo = eng_force(1,&iatom);
      double fxlo = atoms[iatom].fz;

      atoms[iatom].x = xmid + sigma/100.0;
      atoms[iatom].z = zmid;
      double exhi = eng_force(1,&iatom);
      double fxhi = atoms[iatom].fz;

      printf("HOP: x,z = %g,%g; mid = %g,%g,%g; zlo = %g,%g; zhi = %g,%g; xlo = %g,%g; xhi = %g,%g\n",
	     xtmp,ztmp,
	     emid,fxmid,fzmid,ezlo,fzlo,ezhi,fzhi,exlo,fxlo,exhi,fxhi);
      atoms[iatom].x = xtmp;
      atoms[iatom].z = ztmp;
#endif
    }

    dump();
    ntimestep++;

    double foo = eng_force(1,&iatom);
    printf("PRE-RELAX %d %g %g\n",ntimestep,atoms[iatom].fx,atoms[iatom].fz);
    relax(iatom);
    printf("POST-RELAX %d %g %g\n",ntimestep,atoms[iatom].fx,atoms[iatom].fz);
    dump();

    energy = etotal();
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
    if (rsq < cutcountsq) count++;
  }
  return count;
}

/* ----------------------------------------------------------------------
   relax Z position of IATOM while holding all other atoms fixed
   return new Z coord
   return final eng of config and fx,fz on IATOM
   neighbor list for IATOM already exists
------------------------------------------------------------------------- */

double AppSurf::zrelax(int iatom, double *eng, double *fx, double *fz)
{
  double z_initial = atoms[iatom].z;

  // compute bracketed input for dbrent():
  //   lo < mid < hi, eng(mid) < eng(lo), eng(mid) < eng(hi)
  // strategy:
  //   start z at zhi + sigma
  //   reduce z by delta until gradient is up, this is lo
  //   increase z by delta until gradient is down, this is hi
  //   whichever elo/ehi is smaller becomes mid, epsilon away is new lo/hi

  double z,lo,mid,hi,elo,emid,ehi;
  double delta = sigma/2.0;
  double eps = sigma/100.0;
  
  z = zhi + sigma;
  atoms[iatom].z = z;
  elo = eng_force(1,&iatom);
  while (atoms[iatom].fz <= 0.0) {
    z -= delta;
    atoms[iatom].z = z;
    elo = eng_force(1,&iatom);
  }
  lo = z;

  z += delta;
  atoms[iatom].z = z;
  ehi = eng_force(1,&iatom);
  while (atoms[iatom].fz > 0.0) {
    z += delta;
    atoms[iatom].z = z;
    ehi = eng_force(1,&iatom);
  }
  hi = z;

  if (elo < ehi) {
    mid = lo;
    lo = mid - eps;
  } else {
    mid = hi;
    hi = mid + eps;
  }

  // compute optimal znew via dbrent()
  // return eng,fx,fz at znew to caller

#ifdef DEBUG
  double flo,fmid,fhi;
  atoms[iatom].z = lo;
  elo = eng_force(1,&iatom);
  flo = atoms[iatom].fz;
  atoms[iatom].z = mid;
  emid = eng_force(1,&iatom);
  fmid = atoms[iatom].fz;
  atoms[iatom].z = hi;
  ehi = eng_force(1,&iatom);
  fhi = atoms[iatom].fz;
  printf("DB-PRE: %g %g %g; %g %g %g; %g %g %g\n",
  	 lo,elo,flo,mid,emid,fmid,hi,ehi,fhi);
#endif

  double znew;
  double tmp = dbrent(iatom,lo,mid,hi,TOL,&znew);
  atoms[iatom].z = znew;
  *eng = eng_force(1,&iatom);
  *fx = atoms[iatom].fx;
  *fz = atoms[iatom].fz;

#ifdef DEBUG
  printf("DB-POST: %g %g %g\n",znew,eng,fz);
#endif

  atoms[iatom].z = z_initial;
  return znew;
}

/* ----------------------------------------------------------------------
   add an atom to atom list at location IATOM
------------------------------------------------------------------------- */

void AppSurf::add_atom(int iatom, int id, int type, double x, double z)
{
  if (iatom == maxatom) {
    maxatom += DELTA_ATOM;
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

void AppSurf::add_event(int iatom, int style, double rate, double x, double z)
{
  if (nevents == maxevent) {
    maxevent += DELTA_EVENT;
    events = (OneEvent *) 
      memory->srealloc(events,maxevent*sizeof(OneEvent),"surf:events");
    rates = (double *) 
      memory->srealloc(rates,maxevent*sizeof(double),"surf:rates");
  }

  events[nevents].iatom = iatom;
  events[nevents].style = style;
  events[nevents].x = x;
  events[nevents].z = z;
  rates[nevents] = rate;
  nevents++;
}

/* ----------------------------------------------------------------------
   find a diffusive barrier for IATOM hopping to xhop
   return barrier height
   return xnew,znew = saddle point of barrier
   neighbor list for IATOM already exists
------------------------------------------------------------------------- */

double AppSurf::find_barrier(int iatom, double xhop,
			     double *xnew, double *znew)
{
  double eng_initial,x_initial,z_initial;
  double lo,mid,hi;
  double z,eng,fx,fz;

  // compute initial energy

  eng_initial = eng_force(1,&iatom);

  x_initial = atoms[iatom].x;
  z_initial = atoms[iatom].z;

  lo = MIN(x_initial,xhop);
  hi = MAX(x_initial,xhop);

  int niterate = 0;
  while (niterate < 15) {
    mid = 0.5 * (lo+hi);
    atoms[iatom].x = mid;
    z = zrelax(iatom,&eng,&fx,&fz);
    
    if (fx < 0.0) lo = mid;
    else if (fx > 0.0) hi = mid;
    else break;

    niterate++;
  }

  printf("BARRIER: %g %g %g\n",eng-eng_initial,fx,fz);

  *xnew = mid;
  *znew = z;

#ifdef DEBUG
  atoms[iatom].z = z;
  lo = mid - sigma/100.0;
  atoms[iatom].x = lo;
  double elo = eng_force(1,&iatom);
  double fxlo = atoms[iatom].fx;
  double fzlo = atoms[iatom].fz;
  atoms[iatom].x = mid;
  double emid = eng_force(1,&iatom);
  double fxmid = atoms[iatom].fx;
  double fzmid = atoms[iatom].fz;
  hi = mid + sigma/100.0;
  atoms[iatom].x = hi;
  double ehi = eng_force(1,&iatom);
  double fxhi = atoms[iatom].fx;
  double fzhi = atoms[iatom].fz;
  printf("BARRIER: E = %g; lo = %g,%g; mid = %g,%g,%g,%g; hi = %g,%g\n",
	 eng-eng_initial,lo,elo,mid,emid,fxmid,fzmid,hi,ehi);
#endif

  atoms[iatom].x = x_initial;
  atoms[iatom].z = z_initial;

  return eng - eng_initial;
}

/* ----------------------------------------------------------------------
   relax x,z coords of atoms after event occurs where IATOM moved
   set new zhi based on atoms whose coords changed
------------------------------------------------------------------------- */

void AppSurf::relax(int iatom)
{
  neighbor(1,&iatom);
  //sd(1,&iatom);
  cg(1,&iatom);

  zhi = MAX(zhi,atoms[iatom].z);
}

/* ----------------------------------------------------------------------
   return total energy of system
   inefficient double loop over all atoms
------------------------------------------------------------------------- */

double AppSurf::etotal()
{
  int i,j;
  double dx,dz,rsq,r2inv,r6inv;

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
   build neighbor list for N atoms in list within cutoff+skin
   constructed via inefficient double loop
   zero forces on N atoms and atoms in list
------------------------------------------------------------------------- */

void AppSurf::neighbor(int n, int *list)
{
  int i,j,m,jorig;
  int *neighptr;

  if (n > maxatomneigh) {
    delete [] firstneigh;
    delete [] numneigh;
    maxatomneigh = n;
    firstneigh = new int*[maxatomneigh];
    numneigh = new int[maxatomneigh];
  }

  // flag atoms in list
  // only flag owned atoms, ghosts handled in neigh list build

  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++) flag[i] = 0;
  for (i = 0; i < n; i++) flag[list[i]] = 1;

  // build neigh list
  // if atom J is not flagged, always add to list
  // if atom J is flagged, only add if ID of J > ID of I
  // if J is ghost atom, use flag,ID of original atom for testing

  double dx,dz,rsq;

  int npnt = 0;

  for (int ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    atoms[i].fx = atoms[i].fz = 0.0;
    neighptr = &neigh[npnt];
    m = 0;

    for (j = 0; j < nlocal+nghost; j++) {
      if (j == i) continue;
      if (j >= nlocal) jorig = atoms[j].id - 1;
      else jorig = j;
      if (flag[jorig] && atoms[jorig].id < atoms[i].id) continue;

      dx = atoms[i].x - atoms[j].x;
      dz = atoms[i].z - atoms[j].z;
      rsq = dx*dx + dz*dz;
      if (rsq > cutneighsq) continue;

      if (m == maxneigh) {
	maxneigh += DELTA_NEIGH;
	neigh = (int *)
	  memory->srealloc(neigh,maxneigh*sizeof(int),"surf:neigh");
      }
      
      neigh[m++] = j;
      atoms[j].fx = atoms[j].fz = 0.0;
    }

    firstneigh[ilist] = neighptr;
    numneigh[ilist] = m;
    npnt += m;
  }

  delete [] flag;
}

/* ----------------------------------------------------------------------
   build neighbor list for single IATOM of all atoms with dx < cutoff+extra
   zero forces on IATOM and atoms in list
------------------------------------------------------------------------- */

void AppSurf::neighbor_z(int iatom, double extra)
{
  if (maxatomneigh == 0) {
    delete [] firstneigh;
    delete [] numneigh;
    maxatomneigh = 1;
    firstneigh = new int*[maxatomneigh];
    numneigh = new int[maxatomneigh];
  }

  double x = atoms[iatom].x;
  atoms[iatom].fx = atoms[iatom].fz = 0.0;
  double cut = cutoff + extra;
  double dx;

  int m = 0;
  for (int i = 0; i < nlocal+nghost; i++) {
    if (i == iatom) continue;
    dx = fabs(x - atoms[i].x);
    if (dx > cut) continue;

    if (m == maxneigh) {
      maxneigh += DELTA_NEIGH;
      neigh = (int *)
	memory->srealloc(neigh,maxneigh*sizeof(int),"surf:neigh");
    }

    neigh[m++] = i;
    atoms[i].fx = atoms[i].fz = 0.0;
  }
  firstneigh[0] = neigh;
  numneigh[0] = m;
}

/* ----------------------------------------------------------------------
   compute total eng and force on list of N atoms
   neighbor list must have been previously built
   return total eng of N atoms interacting with each other and extra atoms
   set f on N atoms (f on other atoms was zeroed by neigh list build)
------------------------------------------------------------------------- */

double AppSurf::eng_force(int n, int *list)
{
  int i,j,k,nneigh;
  double dx,dz,rsq,r2inv,r6inv,fforce;
  int *neigh;

  // zero forces on N atoms

  double eng = 0.0;
  for (int ii = 0; ii < n; ii++) {
    i = list[ii];
    atoms[i].fx = atoms[i].fz = 0.0;
  }  

  // compute eng and force on N atoms

  for (int ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    neigh = firstneigh[ilist];
    nneigh = numneigh[ilist];
    for (k = 0; k < nneigh; k++) {
      j = neigh[k];
      dx = atoms[i].x - atoms[j].x;
      dz = atoms[i].z - atoms[j].z;
      rsq = dx*dx + dz*dz;
      if (rsq >= cutsq) continue;
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;
      fforce = r6inv * (lj1*r6inv - lj2) * r2inv;
      atoms[i].fx += dx * fforce;
      atoms[i].fz += dz * fforce;
      atoms[j].fx -= dx * fforce;
      atoms[j].fz -= dz * fforce;
      eng += r6inv * (lj3*r6inv - lj4) - offset;
    }
  }

  return eng;
}

/* ----------------------------------------------------------------------
   acquire ghost atoms within cutoff+skin of system boundaries
------------------------------------------------------------------------- */

void AppSurf::ghost_comm()
{
  nghost = 0;
  for (int i = 0; i < nlocal; i++) {
    if (atoms[i].x < xlo+cutneigh) {
      add_atom(nlocal+nghost,atoms[i].id,atoms[i].type,
	       atoms[i].x + xprd,atoms[i].z);
      nghost++;
    } else if (atoms[i].x > xhi-cutneigh) {
      add_atom(nlocal+nghost,atoms[i].id,atoms[i].type,
	       atoms[i].x - xprd,atoms[i].z);
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
  fw = fv = fx = eng_force(1,&iatom);
  dw = dv = dx = -atoms[iatom].fz;

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
      fu = eng_force(1,&iatom);
    } else {
      u = x + SIGN(tol1,d);
      atoms[iatom].z = u;
      fu = eng_force(1,&iatom);
      if (fu > fx) {
	*xmin = x;
	return fx;
      }
    }
    du = -atoms[iatom].fz;

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

/* ----------------------------------------------------------------------
   minimization via steepest descent
------------------------------------------------------------------------- */

void AppSurf::sd(int n, int *list)
{
  int i,ilist,m,fail,niter,neval;
  double alpha,dot;
  double einitial,eprevious,ecurrent;

  int ndof = 2*n;
  double *h = new double[ndof];

  ecurrent = einitial = eng_force(n,list);

  m = 0;
  for (ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    h[m++] = atoms[i].fx;
    h[m++] = atoms[i].fz;
  }

  neval = 0;

  for (niter = 0; niter < MAX_ITER; niter++) {

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    //printf("DIR: %g %g\n",h[0],h[1]);
    //printf("XINIT: eng = %g; %g %g\n",ecurrent,
    //         atoms[list[0]].x,atoms[list[0]].z);
    fail = linemin_scan(n,list,h,ecurrent,DMIN,DMAX,alpha,neval);
    //printf("XFINAL: eng = %g; %g %g\n",ecurrent,
    //	   atoms[list[0]].x,atoms[list[0]].z);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently, all done

    if (neval >= MAX_EVAL) break;

    if (fail || fabs(ecurrent-eprevious) <= 
    	TOLERANCE * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS))
      break;

    // set h to new f = -Grad(x)
    // done if size sq of grad vector < EPS

    dot = 0.0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      dot += atoms[i].fx*atoms[i].fx + atoms[i].fz*atoms[i].fz;
    }
    if (dot < EPS) break;

    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      h[m++] = atoms[i].fx;
      h[m++] = atoms[i].fz;
    }
  }

  delete [] h;

  //printf("SD: %d %d; %g %g\n",niter,neval,einitial,ecurrent);
}

/* ----------------------------------------------------------------------
   minimization of N atoms in list via conjugate gradient iterations
   Polak-Ribiere formulation
------------------------------------------------------------------------- */

void AppSurf::cg(int n, int *list)
{
  int i,ilist,m,gradsearch,fail,niter,neval;
  double alpha,beta,gg,fdotf,fdotg;
  double einitial,eprevious,ecurrent;

  int ndof = 2*n;
  double *g = new double[ndof];
  double *h = new double[ndof];

  ecurrent = einitial = eng_force(n,list);

  m = 0;
  for (ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    h[m] = g[m] = atoms[i].fx;
    m++;
    h[m] = g[m] = atoms[i].fz;
    m++;
  }

  gg = 0.0;
  for (ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    gg += atoms[i].fx*atoms[i].fx + atoms[i].fz*atoms[i].fz;
  }

  neval = 0;
  gradsearch = 1;

  for (niter = 0; niter < MAX_ITER; niter++) {

    // line minimization along direction h from current atom->x

    eprevious = ecurrent;
    fail = linemin_scan(n,list,h,ecurrent,DMIN,DMAX,alpha,neval);

    // if max_eval exceeded, all done
    // if linemin failed or energy did not decrease sufficiently:
    //   all done if searched in grad direction
    //   else force next search to be in grad direction (CG restart)

    if (neval >= MAX_EVAL) break;

    printf("BREAK: %d %g %g %g %g\n",
	   fail,ecurrent,eprevious,ecurrent-eprevious,
	   TOLERANCE * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS));
    if (fail || fabs(ecurrent-eprevious) <= 
    	TOLERANCE * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS)) {
      if (gradsearch == 1) break;
      gradsearch = -1;
    }

    // update h from new f = -Grad(x) and old g
    // old g,h must have migrated with atoms to do this correctly
    // done if size sq of grad vector < EPS
    // force new search dir to be grad dir if need to restart CG
    // set gradsearch to 1 if will search in grad dir on next iteration

    fdotf = fdotg = 0.0;
    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      fdotf += atoms[i].fx*atoms[i].fx + atoms[i].fz*atoms[i].fz;
      fdotg += atoms[i].fx*g[m] + atoms[i].fz*g[m+1];
      m += 2;
    }

    beta = MAX(0.0,(fdotf - fdotg)/gg);
    gg = fdotf;
    if (gg < EPS) break;

    if (gradsearch == -1) beta = 0.0;
    if (beta == 0.0) gradsearch = 1;
    else gradsearch = 0;

    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      g[m] = atoms[i].fx;
      h[m] = g[m] + beta*h[m];
      m++;
      g[m] = atoms[i].fz;
      h[m] = g[m] + beta*h[m];
      m++;
    }
  }

  delete [] g;
  delete [] h;

  printf("CG: %d %d; %g %g\n",niter,neval,einitial,ecurrent);
}

/* ----------------------------------------------------------------------
   line minimization methods
   find minimum-energy starting at x along dir direction
   input: n = # of degrees of freedom on this proc
          x = ptr to atom->x[0] as vector
	  dir = search direction as vector
	  eng = current energy at initial x
	  min/max dist = min/max distance to move any atom coord
   output: return 0 if successful move, set alpha
           return 1 if failed, no move, no need to set alpha
           alpha = distance moved along dir to set x to min-eng config
           caller has several quantities set via last call to eng_force()
	     INSURE last call to eng_force() is consistent with returns
	       if fail, eng_force() of original x
	       if succeed, eng_force() at x + alpha*dir
             atom->x = coords at new configuration
	     atom->f = force (-Grad) is evaulated at new configuration
	     ecurrent = energy of new configuration
   NOTE: when call eng_force: n,x,dir,eng may change due to atom migration
	 updated values are returned by eng_force()
	 this routine CANNOT store atom-based quantities b/c of migration
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   linemin: scan forward by larger and larger steps (SCAN_FACTOR)
   uses no gradient info, but should be very robust
   start at mindist, continue until maxdist
   quit as soon as energy starts to rise
------------------------------------------------------------------------- */

int AppSurf::linemin_scan(int n, int *list, double *dir, double &eng,
			  double mindist, double maxdist,
			  double &alpha, int &nfunc)
{
  int i,ilist,m;
  double fmax,fme,elowest,alphamin,alphamax,alphalast;
  int ndof = 2*n;

  // alphamin = step that moves some atom coord by mindist
  // alphamax = step that moves some atom coord by maxdist

  fme = 0.0;
  for (i = 0; i < n; i++) fme = MAX(fme,fabs(dir[i]));
  //MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  fmax = fme;
  if (fmax == 0.0) return 1;

  alphamin = mindist/fmax;
  alphamax = maxdist/fmax;

  // if minstep is already uphill, fail
  // if eng increases, stop and return previous alpha
  // if alphamax, stop and return alphamax

  elowest = eng;
  alpha = alphamin;

  while (1) {
    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      atoms[i].x += alpha * dir[m++];
      atoms[i].z += alpha * dir[m++];
    }

    eng = eng_force(n,list);
    nfunc++;

    if (alpha == alphamin && eng >= elowest) {
      m = 0;
      for (ilist = 0; ilist < n; ilist++) {
	i = list[ilist];
	atoms[i].x -= alpha * dir[m++];
	atoms[i].z -= alpha * dir[m++];
      }

      eng = eng_force(n,list);
      nfunc++;
      return 1;
    }

    if (eng > elowest) {
      m = 0;
      for (ilist = 0; ilist < n; ilist++) {
	i = list[ilist];
	atoms[i].x += (alphalast-alpha) * dir[m++];
	atoms[i].z += (alphalast-alpha) * dir[m++];
      }

      eng = eng_force(n,list);
      nfunc++;
      alpha = alphalast;
      return 0;
    }      
    if (alpha == alphamax) return 0;

    elowest = eng;
    alphalast = alpha;
    alpha *= SCAN_FACTOR;
    if (alpha > alphamax) alpha = alphamax;
  }
}

/* ----------------------------------------------------------------------
   linemin: use secant approximation to estimate parabola minimum at each step
   should converge more quickly/accurately than "scan", but may be less robust
   how to prevent func evals at new x bigger than maxdist?
   initial secant from two points: 0 and sigma0 = mindist
------------------------------------------------------------------------- */

int AppSurf::linemin_secant(int n, int *list, double *dir, double &eng,
			    double mindist, double maxdist,
			    double &alpha, int &nfunc)
{
  int i,m,ilist,iter;
  double eta,eta_prev,sigma0,alphadelta,fme,fmax,dsq,e0,tmp;
  double epssq = SECANT_EPS * SECANT_EPS;
  int ndof = 2*n;

  // stopping criterion for secant iterations
  // change this

  fme = 0.0;
  for (i = 0; i < ndof; i++) fme += dir[i]*dir[i];
  //MPI_Allreduce(&fme,&dsq,1,MPI_DOUBLE,MPI_SUM,world);
  dsq = fme;

  // sigma0 = smallest allowed step of mindist
  // eval func at sigma0
  // test if minstep is already uphill

  fme = 0.0;
  for (i = 0; i < ndof; i++) fme = MAX(fme,fabs(dir[i]));
  //MPI_Allreduce(&fme,&fmax,1,MPI_DOUBLE,MPI_MAX,world);
  fmax = fme;
  if (fmax == 0.0) return 1;

  sigma0 = mindist/fmax;

  e0 = eng;
  m = 0;
  for (ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    atoms[i].x += sigma0*dir[m++];
    atoms[i].z += sigma0*dir[m++];
  }
  eng = eng_force(n,list);
  printf("  First eng: %g %g\n",sigma0,eng);
  nfunc++;

  if (eng >= e0) {
    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      atoms[i].x -= sigma0*dir[m++];
      atoms[i].z -= sigma0*dir[m++];
    }
    eng = eng_force(n,list);
    nfunc++;
    return 1;
  }

  tmp = 0.0;
  m = 0;
  for (ilist = 0; ilist < n; ilist++) {
    i = list[ilist];
    tmp -= atoms[i].fx * dir[m++];
    tmp -= atoms[i].fz * dir[m++];
  }
  //MPI_Allreduce(&tmp,&eta_prev,1,MPI_DOUBLE,MPI_SUM,world);
  eta_prev = tmp;

  // secant iterations
  // alphadelta = new increment to move, alpha = accumulated move

  alpha = sigma0;
  alphadelta = -sigma0;

  for (iter = 0; iter < MAX_LINE; iter++) {
    alpha += alphadelta;
    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      atoms[i].x += alphadelta * dir[m++];
      atoms[i].z += alphadelta * dir[m++];
    }

    eng = eng_force(n,list);
    printf("  Iter: %g %g\n",alpha,eng);
    nfunc++;

    tmp = 0.0;
    m = 0;
    for (ilist = 0; ilist < n; ilist++) {
      i = list[ilist];
      tmp -= atoms[i].fx * dir[m++];
      tmp -= atoms[i].fz * dir[m++];
    }
    //MPI_Allreduce(&tmp,&eta,1,MPI_DOUBLE,MPI_SUM,world);
    eta = tmp;

    printf("  Alpha adjust: %g %g %g %g\n",
	   alphadelta,eta,eta_prev,eta/(eta_prev-eta));
    alphadelta *= eta / (eta_prev - eta);
    eta_prev = eta;
    if (alphadelta*alphadelta*dsq <= epssq) break;
  }

  // if exited loop on first iteration, func eval was at alpha = 0.0
  // else successful line search

  if (iter == 0) return 1;
  return 0;
}
