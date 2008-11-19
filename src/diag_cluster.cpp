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
#include "diag_cluster.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "random_park.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagCluster::DiagCluster(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  cluster_ids = NULL;
  comm = NULL;
  fp = NULL;
  fpdump = NULL;
  clustlist = NULL;
  opendxroot = NULL;
  ncluster = 0;
  idump = 0;
  dump_style = STANDARD;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fp = fopen(arg[iarg],"w");
	  if (!fp) error->one("Cannot open diag_style cluster output file");
	}
      } else error->all("Illegal diag_style cluster command");
    } else if (strcmp(arg[iarg],"dump") == 0) {
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
		error->one("Cannot open diag_style cluster dump file");
	    }
	  } else error->all("Illegal diag_style cluster command");
	} else if (strcmp(arg[iarg],"opendx") == 0) {
	  idump = 1;
	  dump_style = OPENDX;
	  iarg++;
	  if (iarg < narg) {
	    int n = strlen(arg[iarg]) + 1;
	    opendxroot = new char[n];
	    strcpy(opendxroot,arg[iarg]);
	    opendxcount = 0;
	  } else error->all("Illegal diag_style cluster command");
	} else if (strcmp(arg[iarg],"none") == 0) {
	  idump = 0;
	} else error->all("Illegal diag_style cluster command");
      } else error->all("Illegal diag_style cluster command");
    } else error->all("Illegal diag_style cluster command");
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

DiagCluster::~DiagCluster()
{
  delete comm;
  memory->destroy_1d_T_array(cluster_ids,0);
  free_clustlist();
  delete [] opendxroot;
  if (me == 0 ) {
    if (fp) fclose(fp);
    if (fpdump) fclose(fpdump);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster::init(double time)
{
  if (app->appclass != App::LATTICE)
    error->all("Diag style incompatible with app style");

  applattice = (AppLattice *) app;
  nglobal = applattice->nglobal;
  nlocal = applattice->nlocal;
  nghost = applattice->nghost;
  boxxlo = applattice->boxxlo;
  boxxhi = applattice->boxxhi;
  boxylo = applattice->boxylo;
  boxyhi = applattice->boxyhi;
  boxzlo = applattice->boxzlo;
  boxzhi = applattice->boxzhi;
  xyz = applattice->xyz;
  id = applattice->id;

  memory->destroy_1d_T_array(cluster_ids,0);
  memory->create_1d_T_array(cluster_ids,0,nlocal+nghost-1,"diagcluster:cluster");

  if (!comm) comm = new CommLattice(spk);
  comm->init(NULL,1,0,cluster_ids);

  if (dump_style == OPENDX) {
    if (applattice->latstyle == AppLattice::SC_6N || 
	applattice->latstyle == AppLattice::SC_26N) {
      nx_global = applattice->nx;
      ny_global = applattice->ny;
      nz_global = applattice->nz;
    } else if (applattice->latstyle == AppLattice::SQ_4N || 
	applattice->latstyle == AppLattice::SQ_8N) {
      nx_global = applattice->nx;
      ny_global = applattice->ny;
      nz_global = 1;
    } else {
      error->all("Diag dump_style incompatible with latstyle");
    }
  }

  write_header();
  analyze_clusters(time);

  setup_time(time);
}


/* ---------------------------------------------------------------------- */

void DiagCluster::compute(double time, int iflag, int done)
{
  if (diag_delta > 0.0) iflag = check_time(time, done);

  if (iflag || done) {
    applattice->comm->all();
    analyze_clusters(time);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCluster::analyze_clusters(double time)
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

void DiagCluster::write_header()
{
  if (me == 0) {
    if (fp) {
      fprintf(fp,"Clustering Analysis for Lattice (diag_style cluster) \n");
      fprintf(fp,"nglobal = %d \n",nglobal);
      fprintf(fp,"nprocs = %d \n",nprocs);
    }
  }
}

/* ----------------------------------------------------------------------
   Perform cluster analysis using a definition
   of connectivity provided by the child application
------------------------------------------------------------------------- */

void DiagCluster::generate_clusters()
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
  //     loop m over neighbor pixels of i:
  //       if (spin[m] == spin[n] and not outside domain) push(m)
  //   }

  // void push(m)
  //   id[m] = id
  //   stack[nstack++] = m

  // int pop()
  //   return stack[nstack--]
  //   area++

  // Assign unique global id to each cluster

  // use MPI_Scan to find id_offset
  // loop i over all clusters {
  //    clustlist[i].id += id_offset
  // loop n over all owned sites {
  //    id[n] += id_offset

  // communicate(owned boundary site ids)

  // loop n over all owned boundary sites {
  //     loop m over neighbor pixels of i:
  //       if (spin[m] == spin[n] and m outside domain)
  //          clustlist[id[n]-id_offset].addneigh(id[m])

  // At this point, the problem is reduced to the simpler problem of 
  // clustering the clusters. This can be done by the root process.


  // Update ghost spins
    applattice->comm->all();


  // Set ghost site ids to -1
  for (int i = nlocal; i < nlocal+nghost; i++) {
    cluster_ids[i] = -1;
  }

  // Set local site ids to zero 
  for (int i = 0; i < nlocal; i++) {
    cluster_ids[i] = 0;
  }

  int nclustertot,ii,jj,kk,id;
  double vol,volsum,voltot;

  ncluster = 0;
  volsum = 0.0;

  // loop over all owned sites
  for (int i = 0; i < nlocal; i++) {
    // If already visited, skip
    if (cluster_ids[i] != 0) {
      continue;
    }
    
    // Push first site onto stack
    id = ncluster+1;
    vol = 0.0;
    add_cluster(id,vol,0,NULL);
    cluststack.push(i);
    cluster_ids[i] = id;

    while (cluststack.size()) {
      // First top then pop
      ii = cluststack.top();
      cluststack.pop();
      vol++;
      applattice->push_connected_neighbors(ii,cluster_ids,ncluster,&cluststack);
    }
    clustlist[ncluster-1].volume = vol;
    volsum+=vol;
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
  for (int i = 0; i < nlocal; i++) {
    cluster_ids[i] = clustlist[cluster_ids[i]-1].global_id;
  }

  // Communicate side ids
  comm->all();

  // loop over all owned sites adjacent to boundary
  for (int i = 0; i < nlocal; i++) {
    applattice->connected_ghosts(i,cluster_ids,clustlist,idoffset);
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

void DiagCluster::add_cluster(int id, double vol, int nn, double* neighs)
{
  // grow cluster array

  ncluster++;
  clustlist = (Cluster *) memory->srealloc(clustlist,ncluster*sizeof(Cluster),
					 "diagcluster:clustlist");
  clustlist[ncluster-1] = Cluster(id,vol,nn,neighs);
}

/* ----------------------------------------------------------------------
   dump a snapshot of cluster identities to file pointer fpdump
------------------------------------------------------------------------- */

void DiagCluster::dump_clusters(double time)
{
  int nsend,nrecv;
  double* dbuftmp;
  int maxbuftmp;
  int cid, isite;
  int nsites = nx_global*ny_global*nz_global;
  int* datadx;
  int* randomkeys;
  RandomPark *randomtmp;

  if (me == 0) {
    if (dump_style == STANDARD) {
      fprintf(fpdump,"ITEM: TIMESTEP\n");
      fprintf(fpdump,"%g\n",time);
      fprintf(fpdump,"ITEM: NUMBER OF ATOMS\n");
      fprintf(fpdump,"%d\n",nglobal);
      fprintf(fpdump,"ITEM: BOX BOUNDS\n");
      fprintf(fpdump,"%g %g\n",boxxlo,boxxhi);
      fprintf(fpdump,"%g %g\n",boxylo,boxyhi);
      fprintf(fpdump,"%g %g\n",boxzlo,boxzhi);
      fprintf(fpdump,"ITEM: ATOMS\n");
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
	error->one("Diag style cluster dump file name too long");
      strcpy(filetmp,opendxroot);
      sprintf(filetmp+lroot,"%05d",opendxcount);
      sprintf(filetmp+lroot+lnum,"%s",".dx");
      if (me == 0) {
	fpdump = fopen(filetmp,"w");
	if (!fpdump) error->one("Cannot open diag style cluster dump file");
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
      datadx = (int *) memory->smalloc(nsites*sizeof(int),"diagcluster:datadx");
      randomkeys = (int *) memory->smalloc(ncluster*sizeof(int),"diagcluster:randomkeys");
      randomtmp = new RandomPark(12345);
      for (int i = 0; i < ncluster; i++) {
	randomkeys[i] = randomtmp->irandom(nsites);
      }
    }
  }

  int size_one = 5;

  maxbuftmp = 0;
  MPI_Allreduce(&nlocal,&maxbuftmp,1,MPI_INT,MPI_MAX,world);
  dbuftmp = (double*) memory->smalloc(size_one*maxbuftmp*sizeof(double),
				      "diagcluster:dump_clusters:buftmp");

  int m = 0;

  // pack my lattice values into buffer

  for (int i = 0; i < nlocal; i++) {
    dbuftmp[m++] = id[i];
    dbuftmp[m++] = cluster_ids[i];
    dbuftmp[m++] = xyz[i][0];
    dbuftmp[m++] = xyz[i][1];
    dbuftmp[m++] = xyz[i][2];
  }

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(dbuftmp,size_one*maxbuftmp,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
	nrecv /= size_one;
      } else nrecv = nlocal;
      
      m = 0;

  // print lattice values

      if (dump_style == STANDARD) {
	for (int i = 0; i < nrecv; i++) {
	  cid = static_cast<int>(dbuftmp[m+1])-1;
	  cid = clustlist[cid].global_id;
	  fprintf(fpdump,"%d %d %g %g %g\n",
		  static_cast<int>(dbuftmp[m]),cid,
		  dbuftmp[m+2],dbuftmp[m+3],dbuftmp[m+4]);
	  m += size_one;
	}
      } else if (dump_style == OPENDX) {
	for (int i = 0; i < nrecv; i++) {
	  isite = static_cast<int>(dbuftmp[m]);
	  cid = clustlist[static_cast<int>(dbuftmp[m+1])].global_id;
	  datadx[isite-1] = cid;
	  m += size_one;
	}
      }
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbuftmp,size_one*nlocal,MPI_DOUBLE,0,0,world);
  }

  memory->sfree(dbuftmp);

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

void DiagCluster::free_clustlist()
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

void DiagCluster::stats(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10d",ncluster_reduced);
}

/* ---------------------------------------------------------------------- */

void DiagCluster::stats_header(char *strtmp) {
  if (stats_flag == 0) return;
  sprintf(strtmp," %10s","Nclust");
}
