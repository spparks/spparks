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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_cluster2d.h"
#include "app_lattice2d.h"
#include "comm_lattice2d.h"
#include "random_park.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagCluster2d::DiagCluster2d(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  fp = NULL;
  fpdump = NULL;
  clustlist = NULL;
  opendxroot = NULL;
  ncluster = 0;
  idump = 0;
  dump_style = STANDARD;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style cluster2d output file");
	}
      } else {
	error->all("Illegal diag_style cluster2d command");
      } 
    } else if (strcmp(arg[iarg],"dump_style") == 0) {
      iarg++;
      if (iarg < narg) {
	if (strcmp(arg[iarg],"standard") == 0) {
	  idump = 1;
	  dump_style = STANDARD;
	  iarg++;
	  if (iarg < narg) {
	    if (me == 0) {
	      fpdump = fopen(arg[iarg],"w");
	      if (!fpdump) error->one("Cannot open diag_style cluster2d dump file");
	    }
	  } else {
	    error->all("Illegal diag_style cluster2d command");
	  }
	} else if (strcmp(arg[iarg],"opendx") == 0) {
	  idump = 1;
	  dump_style = OPENDX;
	  iarg++;
	  if (iarg < narg) {
	    int n = strlen(arg[iarg]) + 1;
	    opendxroot = new char[n];
	    strcpy(opendxroot,arg[iarg]);
	    opendxcount = 0;
	  } else {
	    error->all("Illegal diag_style cluster2d command");
	  }
	} else if (strcmp(arg[iarg],"detailed") == 0) {
	  idump = 1;
	  dump_style = DETAILED;
	  iarg++;
	  if (iarg < narg) {
	    if (me == 0) {
	      fpdump = fopen(arg[iarg],"w");
	      if (!fpdump) error->one("Cannot open diag_style cluster2d dump file");
	    }
	  } else {
	    error->all("Illegal diag_style cluster2d command");
	  }
	} else if (strcmp(arg[iarg],"none") == 0) {
	  idump = 0;
	} else {
	    error->all("Illegal diag_style cluster2d command");
	}
      } else {
	error->all("Illegal diag_style cluster2d command");
      }
    } else {
      // error->all("Illegal diag_style cluster2d command");
    }
    iarg++;
  }

}

/* ---------------------------------------------------------------------- */

DiagCluster2d::~DiagCluster2d()
{
  memory->destroy_2d_T_array(cluster_ids,nxlo,nylo);
  free_clustlist();
  delete [] opendxroot;
  if (me == 0 ) {
    if (fp) fclose(fp);
    if (fpdump) fclose(fpdump);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::init(double time)
{
  if (app->appclass != App::LATTICE2D)
    error->all("Diag style incompatible with app style");

  applattice2d = (AppLattice2d *) app;
  nx_global = applattice2d->nx_global;
  ny_global = applattice2d->ny_global;
  nx_procs = applattice2d->nx_procs;
  ny_procs = applattice2d->ny_procs;
  delpropensity = applattice2d->delpropensity;
  nx_local = applattice2d->nx_local;
  ny_local = applattice2d->ny_local;
  nx_offset = applattice2d->nx_offset;
  ny_offset = applattice2d->ny_offset;
  nxlo = applattice2d->nxlo;
  nylo = applattice2d->nylo;
  nxhi = applattice2d->nxhi;
  nyhi = applattice2d->nyhi;

  memory->create_2d_T_array(cluster_ids,nxlo,nxhi,nylo,nyhi,
			    "diagcluster2d:cluster");

  write_header();
  analyze_clusters(time);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagCluster2d::compute(double time, int iflag, int done)
{
  if (diag_delta > 0.0) {
    iflag = check_time(time, done);
  }

  if (iflag) {
    applattice2d->comm->all(applattice2d->lattice);
    analyze_clusters(time);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::analyze_clusters(double time)
{
  if (me == 0) {
    if (fp) {
      fprintf(fp,"\n\n--------------------------------------------------\n");
      fprintf(fp,"Time = %f \n",time);
    }
  }
  free_clustlist();
  generate_clusters();
  if (idump) {
    if (dump_style == STANDARD || dump_style == OPENDX) {
      dump_clusters(time);
    } else if (dump_style == DETAILED) {
      dump_clusters_detailed(time);
    }
  }
}
/* ---------------------------------------------------------------------- */

void DiagCluster2d::write_header()
{
  if (me == 0) {
    if (fp) {
      fprintf(fp,"Clustering Analysis for 2D Lattice (diag_style cluster2d) \n");
      fprintf(fp,"nx_global = %d \n",nx_global);
      fprintf(fp,"ny_global = %d \n",ny_global);
      fprintf(fp,"nprocs = %d \n",nprocs);
    }
  }
}


/* ----------------------------------------------------------------------
   Perform cluster analysis using a definition
   of connectivity provided by the child application
------------------------------------------------------------------------- */

void DiagCluster2d::generate_clusters()
{

  // Psuedocode
  //
  // work on copy of spin array since spin values are changed
  // all id values = 0 initially

  // loop n over all owned sites {
  //   if (id[n] = 0) {
  //     area = 0
  //   } else continue

  //   push(n)
  //   while (stack not empty) {
  //     i = pop()
  //     loop over m neighbor pixels of i:
  //       if (spin[m] == spin[n] and not outside domain) push(m)
  //   }

  // void push(m)
  //   id[m] = id
  //   stack[nstack++] = m

  // int pop()
  //   return stack[nstack--]
  //   area++

  // Set ghost site ids to -1
  // Use four equal sized rectangles arranged in a ring
  for (int j = 1-delpropensity; j <= ny_local; j++) {
    for (int i = 1-delpropensity; i <= 0; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = ny_local+1; j <= ny_local+delpropensity; j++) {
    for (int i = 1-delpropensity; i <= nx_local; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = 1; j <= ny_local+delpropensity; j++) {
    for (int i = nx_local+1; i <= nx_local+delpropensity; i++) {
      cluster_ids[i][j] = -1;
    }
  }
  for (int j = 1-delpropensity; j <= 0; j++) {
    for (int i = 1; i <= nx_local+delpropensity; i++) {
      cluster_ids[i][j] = -1;
    }
  }

  // Set local site ids to zero 
  for (int j = 1; j <= ny_local; j++) {
    for (int i = 1; i <= nx_local; i++) {
      cluster_ids[i][j] = 0;
    }
  }

  int nclustertot,ii,jj,id;
  double vol,volsum,voltot;

  ncluster = 0;
  volsum = 0.0;

  // loop over all owned sites
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {

      // If already visited, skip
      if (cluster_ids[i][j] != 0) {
	continue;
      }

      // Push first site onto stack
      id = ncluster+1;
      vol = 0.0;
      add_cluster(id,vol,0,NULL);
      cluststack.push(i);
      cluststack.push(j);
      cluster_ids[i][j] = id;

      while (cluststack.size()) {
	// First top then pop
	jj = cluststack.top();
	cluststack.pop();
	ii = cluststack.top();
	cluststack.pop();
	vol+=1.0;
	applattice2d->push_connected_neighbors(ii,jj,cluster_ids,ncluster,&cluststack);
      }
      clustlist[ncluster-1].volume = vol;
      volsum+=vol;
    }
  }

  int idoffset;
  MPI_Allreduce(&volsum,&voltot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&ncluster,&nclustertot,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ncluster,&idoffset,1,MPI_INT,MPI_SUM,world);
  idoffset = idoffset-ncluster+1;
  for (int i = 0; i < ncluster; i++) {
    clustlist[i].global_id = i+idoffset;
  }

  // change site ids to global ids
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {
      cluster_ids[i][j] = clustlist[cluster_ids[i][j]-1].global_id;
    }
  }

  applattice2d->comm->all(cluster_ids);

  // loop over all owned sites adjacent to boundary
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {

      applattice2d->connected_ghosts(i,j,cluster_ids,clustlist,idoffset);

    }
  }

  // pack my clusters into buffer
  int me_size,m,maxbuf;
  double* dbufclust;
  int tmp,nrecv;
  MPI_Status status;
  MPI_Request request;
  int nn;
  
  me_size = 0;
  for (int i = 0; i < ncluster; i++) {
    me_size += 3+clustlist[i].nneigh;
  }
  if (me == 0) me_size = 0;

  MPI_Allreduce(&me_size,&maxbuf,1,MPI_INT,MPI_MAX,world);

  dbufclust = new double[maxbuf];

  if (me != 0) {
    m = 0;
    for (int i = 0; i < ncluster; i++) {
      dbufclust[m++] = clustlist[i].global_id;
      dbufclust[m++] = clustlist[i].volume;
      dbufclust[m++] = clustlist[i].nneigh;
      for (int j = 0; j < clustlist[i].nneigh; j++) {
	dbufclust[m++] = clustlist[i].neighlist[j];
      }
    }
    
    if (me_size != m) {
      error->one("Mismatch in counting for dbufclust");
    }

  }

  // proc 0 pings each proc, receives it's data, adds it to list
  // all other procs wait for ping, send their data to proc 0

  if (me == 0) {
    for (int iproc = 1; iproc < nprocs; iproc++) {
      MPI_Irecv(dbufclust,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
      MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
      
      m = 0;
      while (m < nrecv) {
	id = static_cast<int> (dbufclust[m++]);
	vol = dbufclust[m++];
	nn = static_cast<int> (dbufclust[m++]);
	add_cluster(id,vol,nn,&dbufclust[m]);
	m+=nn;
	volsum+=vol;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbufclust,me_size,MPI_DOUBLE,0,0,world);
  }

  delete [] dbufclust;

  // Perform cluster analysis on the clusters

  if (me == 0) {
    int* neighs;
    int jneigh;

    volsum = 0.0;
    ncluster_reduced = 0;

    // loop over all clusters
    for (int i = 0; i < ncluster; i++) {
      
      // If already visited, skip
      if (clustlist[i].volume == 0.0) {
	continue;
      }
      
      // Push first cluster onto stack
      id = clustlist[i].global_id;
      vol = 0.0;
      ncluster_reduced++;
      
      cluststack.push(i);
      vol+=clustlist[i].volume;
      clustlist[i].volume = 0.0;
      
      while (cluststack.size()) {
	// First top then pop
	ii = cluststack.top();
	cluststack.pop();
	
	neighs = clustlist[ii].neighlist;
	for (int j = 0; j < clustlist[ii].nneigh; j++) {
	  jneigh = neighs[j]-idoffset;
	  if (clustlist[jneigh].volume != 0.0) {
	    cluststack.push(jneigh);
	    vol+=clustlist[jneigh].volume;
	    clustlist[jneigh].global_id = id;
	    clustlist[jneigh].volume = 0.0;
	  }
	}
      }
      clustlist[i].volume = vol;
      volsum+=vol;
    }
    
    if (fp) {
      fprintf(fp,"ncluster = %d \nsize = ",ncluster_reduced);
      for (int i = 0; i < ncluster; i++) {
	if (clustlist[i].volume > 0.0) {
	  fprintf(fp," %g",clustlist[i].volume);
	}
      }
      fprintf(fp,"\n");
    }
  }
  

}


/* ---------------------------------------------------------------------- */

void DiagCluster2d::add_cluster(int id, double vol, int nn, double* neighs)
{
  // grow cluster array

  ncluster++;
  clustlist = (Cluster *) memory->srealloc(clustlist,ncluster*sizeof(Cluster),
					 "diagcluster2d:clustlist");
  clustlist[ncluster-1] = Cluster(id,vol,nn,neighs);
}

/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to file pointer fpdump
------------------------------------------------------------------------- */

void DiagCluster2d::dump_clusters(double time)
{
  int nsend,nrecv,nxtmp,nytmp,nxotmp,nyotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;
  int id;
  int isite;
  int nsites = nx_global*ny_global;
  int* datadx;
  int* randomkeys;
  RandomPark *randomtmp;

  if (me == 0) {
    if (dump_style == STANDARD) {
      fprintf(fpdump,"ITEM: TIME\n");
      fprintf(fpdump,"%g\n",time);
      fprintf(fpdump,"ITEM: DIMENSIONS\n");
      fprintf(fpdump,"%d\n%d\n",nx_global,ny_global);
      fprintf(fpdump,"ITEM: ELEMENT CLUSTERID\n");
    } else if (dump_style == OPENDX) {
      int lroot = strlen(opendxroot);
      int lnum;
      int lsuf = 3;
      if (opendxcount == 0) {
	lnum = 1;
      } else {
	lnum = int(log(opendxcount)/log(10))+1;
      }
      if (lnum < 5) lnum = 5;
      char filetmp[100];
      if (99 < lroot+lnum+lsuf)
	error->one("Diag style cluster3d dump file name too long");
      strcpy(filetmp,opendxroot);
      sprintf(filetmp+lroot,"%05d",opendxcount);
      sprintf(filetmp+lroot+lnum,"%s",".dx");
      if (me == 0) {
	fpdump = fopen(filetmp,"w");
	if (!fpdump) error->one("Cannot open diag_style cluster3d dump file");
      }

      opendxcount++;

      fprintf(fpdump,"# Cluster ID dump file for OpenDX\n");
      fprintf(fpdump,"# Time = %g\n",time);
      fprintf(fpdump,"# Create regular grid.\n");
      fprintf(fpdump,"object 1 class gridpositions counts %d %d %d\n",
	      nx_global+1,ny_global+1,2);
      fprintf(fpdump,"origin  0 0 0 \n");
      fprintf(fpdump,"delta   1 0 0 \n");
      fprintf(fpdump,"delta   0 1 0 \n");
      fprintf(fpdump,"delta   0 0 1 \n");
      fprintf(fpdump,"\n# Create connections.\n");
      fprintf(fpdump,"object 2 class gridconnections counts %d %d %d\n",
	      nx_global+1,ny_global+1,2);
      fprintf(fpdump,"\n# Feed data.\n");
      fprintf(fpdump,"object 3 class array type int rank 0 items %d data follows\n#data goes here\n",
	      nx_global*ny_global);
      datadx = (int *) memory->smalloc(nsites*sizeof(int),"diagcluster3d:datadx");
      randomkeys = (int *) memory->smalloc(ncluster*sizeof(int),"diagcluster3d:randomkeys");
      randomtmp = new RandomPark(12345);
      for (int i = 0; i < ncluster; i++) {
	randomkeys[i] = randomtmp->irandom(nsites);
      }
    }
  }
  
  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff
  maxbuftmp = (nx_global/nx_procs+1)*(ny_global/ny_procs+1)+4;
  nsend = nx_local*ny_local+4;
  if (maxbuftmp < nsend) 
    error->one("Maxbuftmp size too small in DiagCluster2d::dump_clusters()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"diagcluster2d:dump_clusters:buftmp");

  // proc 0 writes interactive dump header

  // proc 0 writes timestep header

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;

  // pack my lattice values into buffer
  // Violating normal ordering to satisfy output convention

  for (int j = 1; j <= ny_local; j++) {
    for (int i = 1; i <= nx_local; i++) {
      buftmp[m++] = cluster_ids[i][j];
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
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = buftmp[m++];
      nytmp = buftmp[m++];
      nxotmp = buftmp[m++];
      nyotmp = buftmp[m++];

  // print lattice values
  // isite = global grid cell (1:Nglobal)
  // ordered fast in x, slow in y

      if (dump_style == STANDARD) {
	for (int j = 1; j <= nytmp; j++) {
	  for (int i = 1; i <= nxtmp; i++) {
	    isite = (nyotmp+j-1)*nx_global + (nxotmp+i-1) + 1;
	    id = clustlist[buftmp[m++]-1].global_id;
	    fprintf(fpdump,"%3d %3d \n",isite,id);
	  }
	}
      } else if (dump_style == OPENDX) {
	for (int j = 1; j <= nytmp; j++) {
	  for (int i = 1; i <= nxtmp; i++) {
	    // This ordering of sites matches OpenDX
	    // ordered fast in z, slower in y, slowest in x	      
	    isite = (nxotmp+i-1)*ny_global + nyotmp+j-1;
	    id = randomkeys[clustlist[buftmp[m++]-1].global_id-1];
	    datadx[isite] = id;
	  }
	}
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }
  
  memory->sfree(buftmp);

  if (me == 0) {
    if (dump_style == STANDARD) {
      //
    } else if (dump_style == OPENDX) {
      isite = 0;
      while (isite < nsites) {
	for (int i = 0; i < 20; i++) {
	  fprintf(fpdump,"%d ",datadx[isite++]);
	  if (isite == nsites) break;
	}
	fprintf(fpdump,"\n");
      }
      fprintf(fpdump,"attribute \"dep\" string \"connections\"\n");
      fprintf(fpdump,"\n# Create a field.\n");
      fprintf(fpdump,"object \"9 grain microstructure\" class field\n");
      fprintf(fpdump,"component \"positions\" value 1\n");
      fprintf(fpdump,"component \"connections\" value 2\n");
      fprintf(fpdump,"component \"data\" value 3\n");
      fprintf(fpdump,"\nend\n");
      
      fclose(fpdump);
      fpdump = NULL;
      
      memory->sfree(datadx);
      memory->sfree(randomkeys);
      delete randomtmp;
      
    }
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::free_clustlist()
{
  // Can not call Cluster destructor, because 
  // that would free memory twice.
  // Instead, need to delete neighlist manually.
  for (int i = 0; i < ncluster; i++) {
    free(clustlist[i].neighlist);
  }
  memory->sfree(clustlist);
  clustlist = NULL;
  ncluster = 0;
}

/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to the screen in 2D layout
   all the cluster identities for each processor domain are printed out
------------------------------------------------------------------------- */

void DiagCluster2d::dump_clusters_detailed(double time)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxhtmp,nyhtmp,nzhtmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;
  int id;

  // set up communication buffer
  // maxbuftmp must equal the maximum number of spins on one domain 
  // plus some extra stuff

  maxbuftmp = ((nx_global-1)/nx_procs+1+2*delpropensity)*
    ((ny_global-1)/ny_procs+1+2*delpropensity)+9;
  nsend = (nx_local+2*delpropensity)*(ny_local+2*delpropensity)+9;
  if (maxbuftmp < nsend) 
    error->one("Maxbuftmp size too small in DiagCluster2d::dump_clusters()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),
				  "diagcluster2d:dump_clusters:buftmp");

  // proc 0 writes interactive dump header

  if (me == 0) {
    fprintf(fpdump,"*** Cluster Dump ***\n");
    fprintf(fpdump,"Time = %f \n",time);
    fprintf(fpdump,"nx_global = %d ny_global = %d\n",nx_global,ny_global);
  }

  int m = 0;

  // pack local layout info into buffer

  buftmp[m++] = nx_local;
  buftmp[m++] = ny_local;
  buftmp[m++] = 0;
  // Need to delete these two
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = 0;
  buftmp[m++] = nx_offset;
  buftmp[m++] = ny_offset;
  buftmp[m++] = 0;

  // pack my lattice values into buffer
  // Need to violate normal ordering in order to simplify output

  for (int j = 1-delpropensity; j <= ny_local+delpropensity; j++) {
    for (int i = 1-delpropensity; i <= nx_local+delpropensity; i++) {
      buftmp[m++] = cluster_ids[i][j];
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
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = buftmp[m++];
      nytmp = buftmp[m++];
      nztmp = buftmp[m++];
      nxhtmp = buftmp[m++];
      nyhtmp = buftmp[m++];
      nzhtmp = buftmp[m++];
      nxotmp = buftmp[m++];
      nyotmp = buftmp[m++];
      nzotmp = buftmp[m++];
      fprintf(fpdump,"iproc = %d \n",iproc);
      fprintf(fpdump,"nxlocal = %d \nnylocal = %d \n",nxtmp,nytmp);
      fprintf(fpdump,"nx_offset = %d \nny_offset = %d \n",nxotmp,nyotmp);
      m = nrecv;
      for (int j = nytmp+delpropensity; j >= 1-delpropensity; j--) {
	m-=nxtmp+2*delpropensity;
	for (int i = 1-delpropensity; i <= nxtmp+delpropensity; i++) {
	  // work out the correct global id of this site
	  id = clustlist[buftmp[m++]-1].global_id;
	  fprintf(fpdump,"%3d",id);
	}
	fprintf(fpdump,"\n");
	m-=nxtmp+2*delpropensity;
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buftmp,nsend,MPI_INT,0,0,world);
  }

  memory->sfree(buftmp);
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10d",ncluster_reduced);
}

/* ---------------------------------------------------------------------- */

void DiagCluster2d::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Nclust");
}
