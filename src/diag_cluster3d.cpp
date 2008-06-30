/* -----------------------q-----------------------------------------------
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
#include "diag_cluster3d.h"
#include "app_lattice3d.h"
#include "comm_lattice3d.h"
#include "random_park.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagCluster3d::DiagCluster3d(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  fp = NULL;
  fpdump = NULL;
  clustlist = NULL;
  opendxroot = NULL;
  ncluster = 0;
  idump = 0;
  dump_style = STANDARD;
  radius = 0.0;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style cluster3d output file");
	}
      } else {
	error->all("Illegal diag_style cluster3d command");
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
	      if (!fpdump)
		error->one("Cannot open diag_style cluster3d dump file");
	    }
	  } else {
	    error->all("Illegal diag_style cluster3d command");
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
	    error->all("Illegal diag_style cluster3d command");
	  }
	} else if (strcmp(arg[iarg],"none") == 0) {
	  idump = 0;
	} else {
	    error->all("Illegal diag_style cluster3d command");
	}
      } else {
	error->all("Illegal diag_style cluster3d command");
      }
    } else {
      //      error->all("Illegal diag_style cluster3d command");
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagCluster3d::~DiagCluster3d()
{
  memory->destroy_3d_T_array(cluster_ids,nxlo,nylo,nzlo);
  free_clustlist();
  delete [] opendxroot;
  if (me == 0 ) {
    if (fp) fclose(fp);
    if (fpdump) fclose(fpdump);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster3d::init(double time)
{
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

  memory->create_3d_T_array(cluster_ids,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
			    "diagcluster3d:cluster");

  write_header();
  analyze_clusters(time);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagCluster3d::compute(double time, int done)
{
  if (check_time(time, done)) {
    applattice3d->comm->all(applattice3d->lattice);
    analyze_clusters(time);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster3d::analyze_clusters(double time)
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
    dump_clusters(time);
  }
}
/* ---------------------------------------------------------------------- */

void DiagCluster3d::write_header()
{
  if (me == 0) {
    if (fp) {
      fprintf(fp,"Clustering Analysis for 3D Lattice (diag_style cluster3d) \n");
      fprintf(fp,"nx_global = %d \n",nx_global);
      fprintf(fp,"ny_global = %d \n",ny_global);
      fprintf(fp,"nz_global = %d \n",nz_global);
      fprintf(fp,"nprocs = %d \n",nprocs);
    }
  }
}

/* ----------------------------------------------------------------------
   Perform cluster analysis using a definition
   of connectivity provided by the child application
------------------------------------------------------------------------- */

void DiagCluster3d::generate_clusters()
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
  // Use six slabs, with z-slabs grabbing
  // the xz and yz ghost edges and y-slab grabbing
  // the xy edges
  for (int k = 1-delghost; k <= 0; k++) {
    for (int j = 1-delghost; j <= ny_local+delghost; j++) {
      for (int i = 1-delghost; i <= nx_local+delghost; i++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }
  for (int k = nz_local+1; k <= nz_local+delghost; k++) {
    for (int j = 1-delghost; j <= ny_local+delghost; j++) {
      for (int i = 1-delghost; i <= nx_local+delghost; i++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }
  for (int j = 1-delghost; j <= 0; j++) {
    for (int k = 1; k <= nz_local; k++) {
      for (int i = 1-delghost; i <= nx_local+delghost; i++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }
  for (int j = ny_local+1; j <= ny_local+delghost; j++) {
    for (int k = 1; k <= nz_local; k++) {
      for (int i = 1-delghost; i <= nx_local+delghost; i++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }
  for (int i = 1-delghost; i <= 0; i++) {
    for (int k = 1; k <= nz_local; k++) {
      for (int j = 1; j <= ny_local; j++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }
  for (int i = nx_local+1; i <= nx_local+delghost; i++) {
    for (int k = 1; k <= nz_local; k++) {
      for (int j = 1; j <= ny_local; j++) {
	cluster_ids[i][j][k] = -1;
      }
    }
  }

  // Set local site ids to zero 
  for (int k = 1; k <= nz_local; k++) {
    for (int j = 1; j <= ny_local; j++) {
      for (int i = 1; i <= nx_local; i++) {
	cluster_ids[i][j][k] = 0;
      }
    }
  }

  double vol,volsum,voltot;
  int nclustertot,ii,jj,kk,id;

  ncluster = 0;
  volsum = 0.0;

  // loop over all owned sites
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {
      for (int k = 1; k <= nz_local; k++) {

      // If already visited, skip
      if (cluster_ids[i][j][k] != 0) {
	continue;
      }

      // Push first site onto stack
      id = ncluster+1;
      vol = 0.0;
      add_cluster(id,vol,0,NULL);
      cluststack.push(i);
      cluststack.push(j);
      cluststack.push(k);
      cluster_ids[i][j][k] = id;

      while (cluststack.size()) {
	// First top then pop
	kk = cluststack.top();
	cluststack.pop();
	jj = cluststack.top();
	cluststack.pop();
	ii = cluststack.top();
	cluststack.pop();
	vol+=1.0;
	applattice3d->push_connected_neighbors(ii,jj,kk,cluster_ids,ncluster,&cluststack);
      }
      clustlist[ncluster-1].volume = vol;
      volsum+=vol;
      }
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
      for (int k = 1; k <= nz_local; k++) {
	cluster_ids[i][j][k] = clustlist[cluster_ids[i][j][k]-1].global_id;
      }
    }
  }

  applattice3d->comm->all(cluster_ids);

  // loop over all owned sites adjacent to boundary
  for (int i = 1; i <= nx_local; i++) {
    for (int j = 1; j <= ny_local; j++) {
      for (int k = 1; k <= nz_local; k++) {

      applattice3d->connected_ghosts(i,j,k,cluster_ids,clustlist,idoffset);

      }
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

    ncluster_reduced = 0;

    // loop over all clusters
    for (int i = 0; i < ncluster; i++) {
      
      // If already visited, skip
      if (clustlist[i].volume < 0.0) {
	continue;
      }
      
      // Push first cluster onto stack
      id = clustlist[i].global_id;
      vol = 0.0;
      ncluster_reduced++;
      //      fprintf(fp,"vol new = %g \n",clustlist[i].volume);
      cluststack.push(i);
      vol+=clustlist[i].volume;
      clustlist[i].volume = -1.0;
      
      while (cluststack.size()) {
	// First top then pop
	ii = cluststack.top();
	cluststack.pop();
	
	neighs = clustlist[ii].neighlist;
	for (int j = 0; j < clustlist[ii].nneigh; j++) {
	  jneigh = neighs[j]-idoffset;
	  if (clustlist[jneigh].volume > 0.0) {
	    cluststack.push(jneigh);
	    vol+=clustlist[jneigh].volume;
	    clustlist[jneigh].global_id = id;
	    clustlist[jneigh].volume = -1.0;
	  }
	}
      }
      clustlist[i].volume = vol;
    }

    radius = 0.0;
    double onethird = 1.0/3.0;
    for (int i = 0; i < ncluster; i++) {
      if (clustlist[i].volume > 0.0) {
	radius += pow(clustlist[i].volume,onethird);
      }
    }
    radius /= ncluster_reduced;

    if (fp) {
      fprintf(fp,"ncluster_reduced = %d \nvolumes = ",ncluster_reduced);
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

void DiagCluster3d::add_cluster(int id, double vol, int nn, double* neighs)
{
  // grow cluster array

  ncluster++;
  clustlist = (Cluster *) memory->srealloc(clustlist,ncluster*sizeof(Cluster),
					 "diagcluster3d:clustlist");
  clustlist[ncluster-1] = Cluster(id,vol,nn,neighs);
}

/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to file pointer fpdump
------------------------------------------------------------------------- */

void DiagCluster3d::dump_clusters(double time)
{
  int nsend,nrecv,nxtmp,nytmp,nztmp,nxotmp,nyotmp,nzotmp;
  int size_one = 1;
  int* buftmp;
  int maxbuftmp;
  int id;
  int isite;
  int nsites = nx_global*ny_global*nz_global;
  int* datadx;
  int* randomkeys;
  RandomPark *randomtmp;

  if (me == 0) {
    if (dump_style == STANDARD) {
      fprintf(fpdump,"ITEM: TIME\n");
      fprintf(fpdump,"%g\n",time);
      fprintf(fpdump,"ITEM: DIMENSIONS\n");
      fprintf(fpdump,"%d\n%d\n%d\n",nx_global,ny_global,nz_global);
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
	if (!fpdump) error->one("Cannot open diag style cluster3d dump file");
      }

      opendxcount++;

      fprintf(fpdump,"# Cluster ID dump file for OpenDX\n");
      fprintf(fpdump,"# Time = %g\n",time);
      fprintf(fpdump,"# Create regular grid.\n");
      fprintf(fpdump,"object 1 class gridpositions counts %d %d %d\n",
	      nx_global+1,ny_global+1,nz_global+1);
      fprintf(fpdump,"origin  0 0 0 \n");
      fprintf(fpdump,"delta   1 0 0 \n");
      fprintf(fpdump,"delta   0 1 0 \n");
      fprintf(fpdump,"delta   0 0 1 \n");
      fprintf(fpdump,"\n# Create connections.\n");
      fprintf(fpdump,"object 2 class gridconnections counts %d %d %d\n",
	      nx_global+1,ny_global+1,nz_global+1);
      fprintf(fpdump,"\n# Feed data.\n");
      fprintf(fpdump,"object 3 class array type int rank 0 items %d data follows\n#data goes here\n",
	      nx_global*ny_global*nz_global);
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
  maxbuftmp = (nx_global/nx_procs+1)*(ny_global/ny_procs+1)*
    (nz_global/nz_procs+1)+6;
  nsend = nx_local*ny_local*nz_local+6;
  if (maxbuftmp < nsend) 
    error->one("Maxbuftmp size too small in DiagCluster3d::dump_clusters()");
  
  buftmp = (int*) memory->smalloc(maxbuftmp*sizeof(int),"diagcluster3d:dump_clusters:buftmp");

  // proc 0 writes interactive dump header

  // proc 0 writes timestep header

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
	buftmp[m++] = cluster_ids[i][j][k];
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
	MPI_Irecv(buftmp,maxbuftmp,MPI_INT,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_INT,&nrecv);
      } else nrecv = nsend;
      
      m = 0;
      nxtmp = buftmp[m++];
      nytmp = buftmp[m++];
      nztmp = buftmp[m++];
      nxotmp = buftmp[m++];
      nyotmp = buftmp[m++];
      nzotmp = buftmp[m++];

  // print lattice values
  // isite = global grid cell (1:Nglobal)
  // ordered fast in x, slower in y, slowest in z

      if (dump_style == STANDARD) {
	for (int k = 1; k <= nztmp; k++) {
	  for (int j = 1; j <= nytmp; j++) {
	    for (int i = 1; i <= nxtmp; i++) {
	      isite = ((nzotmp+k-1)*ny_global + nyotmp+j-1)
		*nx_global + (nxotmp+i-1) + 1;
	      id = clustlist[buftmp[m++]-1].global_id;
	      fprintf(fpdump,"%3d %3d \n",isite,id);
	    }
	  }
	}
      } else if (dump_style == OPENDX) {
	for (int k = 1; k <= nztmp; k++) {
	  for (int j = 1; j <= nytmp; j++) {
	    for (int i = 1; i <= nxtmp; i++) {
	      // This ordering of sites matches OpenDX
	      // ordered fast in z, slower in y, slowest in x	      
	      isite = ((nxotmp+i-1)*ny_global + nyotmp+j-1)
		*nz_global + (nzotmp+k-1);
	      id = randomkeys[clustlist[buftmp[m++]-1].global_id-1];
	      datadx[isite] = id;
	    }
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

void DiagCluster3d::free_clustlist()
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

/* ---------------------------------------------------------------------- */

void DiagCluster3d::stats(char *strtmp) {
  sprintf(strtmp," %10d %10g",ncluster_reduced,radius);
}

/* ---------------------------------------------------------------------- */

void DiagCluster3d::stats_header(char *strtmp) {
  sprintf(strtmp," %10s %10s","Nclust","AveRadius");
}
