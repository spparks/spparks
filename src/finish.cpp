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

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Finish::Finish (SPPARKS *spk, int flag) : Pointers(spk)
{
  double time,tmp;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // deduce time_other

  double time_other = timer->array[TIME_LOOP] -
    (timer->array[TIME_SOLVE] + timer->array[TIME_COMM] + 
     timer->array[TIME_UPDATE] + timer->array[TIME_OUTPUT] + 
     timer->array[TIME_APP]);

  double time_loop = timer->array[TIME_LOOP];
  MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_loop = tmp/nprocs;

  // overall loop time

  if (me == 0) {
    if (screen) 
      fprintf(screen,
	      "Loop time of %g on %d procs\n",time_loop,nprocs);
    if (logfile)
      fprintf(logfile,
	      "Loop time of %g on %d procs\n",time_loop,nprocs);
  }

  if (flag == 0) return;

  if (me == 0) {
    if (screen) fprintf(screen,"\n");
    if (logfile) fprintf(logfile,"\n");
  }

  // timing breakdowns

  if (time_loop == 0.0) time_loop = 1.0;

  time = timer->array[TIME_SOLVE];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"Solve time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"Solve time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_UPDATE];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"Update time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"Update time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_COMM];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"Comm  time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"Comm  time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_OUTPUT];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"Outpt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"Outpt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_APP];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"App   time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"App   time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = time_other;
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    if (screen) 
      fprintf(screen,"Other time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    if (logfile) 
      fprintf(logfile,"Other time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }
}

/* ---------------------------------------------------------------------- */

void Finish::stats(int n, double *data, 
		   double *pave, double *pmax, double *pmin,
		   int nhisto, int *histo)
{
  int i,m;
  int *histotmp;

  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;
  for (i = 0; i < n; i++) {
    ave += data[i];
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  int ntotal;
  MPI_Allreduce(&n,&ntotal,1,MPI_INT,MPI_SUM,world);
  double tmp;
  MPI_Allreduce(&ave,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  ave = tmp/ntotal;
  MPI_Allreduce(&min,&tmp,1,MPI_DOUBLE,MPI_MIN,world);
  min = tmp;
  MPI_Allreduce(&max,&tmp,1,MPI_DOUBLE,MPI_MAX,world);
  max = tmp;

  for (i = 0; i < nhisto; i++) histo[i] = 0;

  double del = max - min;
  for (i = 0; i < n; i++) {
    if (del == 0.0) m = 0;
    else m = static_cast<int> ((data[i]-min)/del * nhisto);
    if (m > nhisto-1) m = nhisto-1;
    histo[m]++;
  }

  histotmp = (int *) memory->smalloc(nhisto*sizeof(int),
					  "finish:histotmp");
  MPI_Allreduce(histo,histotmp,nhisto,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
  memory->sfree(histotmp);

  *pave = ave;
  *pmax = max;
  *pmin = min;
}
