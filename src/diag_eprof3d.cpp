/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_eprof3d.h"
#include "app_lattice3d.h"
#include "comm_lattice3d.h"

using namespace SPPARKS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEprof3d::DiagEprof3d(SPK *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  fp = NULL;
  prof = NULL;
  ndata = 0;
  prof_style = STANDARD;
  prof_index = 0;
  iboundary = 0;
  eav = 0.0;
  eb1 = 0.0;
  eb2 = 0.0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style eprof3d output file");
	}
      } else {
	error->all("Illegal diag_style cluster3d command");
      } 
    } else if (strcmp(arg[iarg],"axis") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"x") == 0) {
	  prof_index = 0;
	} else if (strcmp(arg[iarg],"y") == 0) {
	  prof_index = 1;
	} else if (strcmp(arg[iarg],"z") == 0) {
	  prof_index = 2;
	} else {
	  error->all("Illegal diag_style eprof3d command");
	}
      } else {
	error->all("Illegal diag_style eprof3d command");
      }
    } else if (strcmp(arg[iarg],"boundary") == 0) {
      iboundary = 1;
    } else {
      //      error->all("Illegal diag_style eprof3d command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagEprof3d::~DiagEprof3d()
{
  memory->sfree(prof);
  memory->sfree(count);
  if (iboundary == 1) {
    memory->sfree(ixtable);
    memory->sfree(iytable);
    memory->sfree(iztable);
  }
  if (me == 0 ) {
    if (fp) fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

void DiagEprof3d::init(double time)
{
  int ntmp;

  if (app->appclass != App::LATTICE3D) error->all("diag_style incompatible with app_style");

  applattice3d = (AppLattice3d *) app;
  nx_global = applattice3d->nx_global;
  ny_global = applattice3d->ny_global;
  nz_global = applattice3d->nz_global;
  nx_procs = applattice3d->nx_procs;
  ny_procs = applattice3d->ny_procs;
  nz_procs = applattice3d->nz_procs;
  delghost = applattice3d->delghost;
  nx_local = applattice3d->nx_local;
  ny_local = applattice3d->ny_local;
  nz_local = applattice3d->nz_local;
  nx_offset = applattice3d->nx_offset;
  ny_offset = applattice3d->ny_offset;
  nz_offset = applattice3d->nz_offset;
  nxlo = applattice3d->nxlo;
  nylo = applattice3d->nylo;
  nzlo = applattice3d->nzlo;
  nxhi = applattice3d->nxhi;
  nyhi = applattice3d->nyhi;
  nzhi = applattice3d->nzhi;

  if (iboundary == 0) {
    if (prof_index == 0) {
      ndata = nx_global;
    } else if (prof_index == 1) {
      ndata = ny_global;
    } else if (prof_index == 2) {
      ndata = nz_global;
    }
  } else {
    nbound = nx_global/(4*nx_procs)+1;
    ntmp = ny_global/(4*ny_procs)+1;
    nbound = MAX(nbound,ntmp);
    ntmp = nz_global/(4*nz_procs)+1;
    nbound = MAX(nbound,ntmp);
    ndata = 2*nbound+1;

    ntmp = nx_global/nx_procs+1;
    ixtable = (int*) memory->smalloc(ntmp*sizeof(int),"diageprof3d:write_prof:ixtable");
    ntmp = ny_global/ny_procs+1;
    iytable = (int*) memory->smalloc(ntmp*sizeof(int),"diageprof3d:write_prof:iytable");
    ntmp = nz_global/nz_procs+1;
    iztable = (int*) memory->smalloc(ntmp*sizeof(int),"diageprof3d:write_prof:iztable");

  }

  prof = (double *) memory->smalloc(ndata*sizeof(double),"diageprof3d:prof");
  count = (double *) memory->smalloc(ndata*sizeof(double),"diageprof3d:count");

  write_header();
  write_prof(time);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagEprof3d::compute(double time, int done)
{
  if (check_time(time, done)) write_prof(time);
}

/* ----------------------------------------------------------------------
   write header for energy profile for snapshot to file pointer fp
------------------------------------------------------------------------- */

void DiagEprof3d::write_header()
{
}

/* ----------------------------------------------------------------------
   write energy profile for snapshot to file pointer fp
------------------------------------------------------------------------- */

void DiagEprof3d::write_prof(double time)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  double* buftmp;
  int maxbuftmp;
  int iprof,iprofz,iprofyz;
  int nxhtmp,nyhtmp,nzhtmp,ix,iy,iz;
  int absx,absy,absz,absyz;

  applattice3d->comm->all(applattice3d->lattice);

  if (me == 0) {
    if (prof_style == STANDARD) {
      fprintf(fp,"ITEM: TIME\n");
      fprintf(fp,"%g\n",time);
      fprintf(fp,"ITEM: NDATA\n");
      fprintf(fp,"%d\n",nx_global);
      fprintf(fp,"ITEM: INDEX ENERGY/SITE\n");
      
      for (int iprof = 0; iprof < ndata; iprof++) {
	prof[iprof] = 0.0;
	count[iprof] = 0.0;
      }

    }
  }

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff
  maxbuftmp = (nx_global/nx_procs+1)*(ny_global/ny_procs+1)*
    (nz_global/nz_procs+1)+6;
  nsend = nx_local*ny_local*nz_local+6;
  if (maxbuftmp < nsend) 
    error->one("maxbuftmp size too small in DiagEprof3d::write_prof()");
  
  buftmp = (double*) memory->smalloc(maxbuftmp*sizeof(double),"diageprof3d:write_prof:buftmp");

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = nz_local;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = nz_offset;

  // pack my lattice values into buffer
  // Violating normal ordering to satisfy output convention

  for (int k = 1; k <= nz_local; k++) {
    for (int j = 1; j <= ny_local; j++) {
      for (int i = 1; i <= nx_local; i++) {
	buftmp[m++] = applattice3d->site_energy(i,j,k);;
      }
    }
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buftmp,maxbuftmp,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = (int) buftmp[m++];
      nytmp = (int) buftmp[m++];
      nztmp = (int) buftmp[m++];
      nxotmp = (int) buftmp[m++];
      nyotmp = (int) buftmp[m++];
      nzotmp = (int) buftmp[m++];

  // sum lattice energies to profile
  // isite = global grid cell (1:Nglobal)
  // ordered fast in x, slower in y, slowest in z

      if (prof_style == STANDARD) {
	if (iboundary == 0) {
	  for (int k = 1; k <= nztmp; k++) {
	    for (int j = 1; j <= nytmp; j++) {
	      for (int i = 1; i <= nxtmp; i++) {
		if (prof_index == 0) {
		  iprof = i+nxotmp-1;
		} else if (prof_index == 1) {
		  iprof = j+nyotmp-1;
		} else if (prof_index == 2) {
		  iprof = k+nzotmp-1;
		}
		prof[iprof]+=buftmp[m++];
		count[iprof]++;
	      }
	    }
	  }
	} else {

	  nxhtmp = nxtmp/2+1;
	  for (int i = 1; i <= nxtmp; i++) {
	    if (i < nxhtmp) {
	      ix = nbound-1-MIN(i-1,nxhtmp-1-i);
	    } else {
	      ix = nbound+1+MIN(i-nxhtmp,nxtmp-i);
	    }
	    ixtable[i] = ix;
	  }

	  nyhtmp = nytmp/2+1;
	  for (int j = 1; j <= nytmp; j++) {
	    if (j < nyhtmp) {
	      iy = nbound-1-MIN(j-1,nyhtmp-1-j);
	    } else {
	      iy = nbound+1+MIN(j-nyhtmp,nytmp-j);
	    }
	    iytable[j] = iy;
	  }

	  nzhtmp = nztmp/2+1;
	  for (int k = 1; k <= nztmp; k++) {
	    if (k < nzhtmp) {
	      iz = nbound-1-MIN(k-1,nzhtmp-1-k);
	    } else {
	      iz = nbound+1+MIN(k-nzhtmp,nztmp-k);
	    }
	    iztable[k] = iz;
	  }

	  for (int k = 1; k <= nztmp; k++) {

	    iprofz = iztable[k]; 
	    absz = abs(iprofz-nbound);

	    for (int j = 1; j <= nytmp; j++) {

	      absy = abs(iytable[j]-nbound);
	      if (absy < absz) {
		iprofyz = iytable[j];
		absyz = absy;
	      } else {
		iprofyz = iprofz;
		absyz = absz;
	      }

	      for (int i = 1; i <= nxtmp; i++) {

		absx = abs(ixtable[i]-nbound);
		if (absx < absyz) {
		  iprof = ixtable[i];
		} else {
		  iprof = iprofyz;
		}

		prof[iprof]+=buftmp[m++];
		count[iprof]++;
	      }
	    }
	  }
	}
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_DOUBLE,0,0,world);
  }

  memory->sfree(buftmp);
  
  if (me == 0) {
    int ii;
    double etmp;

    eav = 0.0;

    if (prof_style == STANDARD) {
      for (int iiprof = 0; iiprof < ndata; iiprof++) {
	if (iboundary == 0) {
	  ii = iiprof+1;
	} else {
	  if (iiprof < nbound) {
	    ii = iiprof - nbound;
	  } else {
	    ii = iiprof - nbound;
	  }
	}

	if (count[iiprof] > 0) {
	  etmp = prof[iiprof]/count[iiprof]	  ;
	  eav += prof[iiprof];
	} else {
	  etmp = 0.0;
	}
	fprintf(fp,"%d %g %g \n",ii,etmp,count[iiprof]);
      }
      eav /= nx_global*ny_global*nz_global;
      iprof = nbound-1;
      if (count[iprof] > 0) eb1 = prof[iprof]/count[iprof];
      else eb1 = 0.0;
      iprof = nbound+1;
      if (count[iprof] > 0) eb2 = prof[iprof]/count[iprof];
      else eb2 = 0.0;
    }
  }

}

/* ---------------------------------------------------------------------- */

void DiagEprof3d::stats(char *strtmp) {
  if (iboundary == 0) {
    sprintf(strtmp," %10g",eav);
  } else {
    sprintf(strtmp," %10g %10g %10g",eav,eb1,eb2);
  }
}

/* ---------------------------------------------------------------------- */

void DiagEprof3d::stats_header(char *strtmp) {
  if (iboundary == 0) {
    sprintf(strtmp," %10s","Eav");
  } else {
    sprintf(strtmp," %10s %10s %10s","Eav","Eb1","Eb2");
  }
}
