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

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "irregular.h"
#include "memory.h"

using namespace SPPARKS_NS;

// allocate space for static class variable
// prototype for non-class function

int *Irregular::proc_recv_copy;
int compare_standalone(const void *, const void *);

enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ---------------------------------------------------------------------- */

Irregular::Irregular(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // migrate work vectors

  maxlocal = 0;
  mproclist = NULL;
  msizes = NULL;

  // send buffers

  maxdbuf = 0;
  dbuf = NULL;
  maxbuf = 0;
  buf = NULL;

  // universal work vectors

  memory->create(work1,nprocs,"irregular:work1");
  memory->create(work2,nprocs,"irregular:work2");

  // initialize buffers for migrate atoms, not used for datum comm
  // these can persist for multiple irregular operations

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ---------------------------------------------------------------------- */

Irregular::~Irregular()
{
  memory->destroy(mproclist);
  memory->destroy(msizes);
  memory->destroy(dbuf);
  memory->destroy(buf);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
   create communication plan based on list of datums of uniform size
   n = # of datums to send
   proclist = proc to send each datum to, can include self
   sortflag = flag for sorting order of received messages by proc ID
   return total # of datums I will recv, including any to self
------------------------------------------------------------------------- */

int Irregular::create_data(int n, int *proclist, int sortflag)
{
  int i,m;

  // setup for collective comm
  // work1 = 1 for procs I send a message to, not including self
  // work2 = 1 for all procs, used for ReduceScatter

  for (i = 0; i < nprocs; i++) {
    work1[i] = 0;
    work2[i] = 1;
  }
  for (i = 0; i < n; i++) work1[proclist[i]] = 1;
  work1[me] = 0;

  // nrecv_proc = # of procs I receive messages from, not including self

  MPI_Reduce_scatter(work1,&nrecv_proc,work2,MPI_INT,MPI_SUM,world);

  // allocate receive arrays

  proc_recv = new int[nrecv_proc];
  num_recv = new int[nrecv_proc];
  request = new MPI_Request[nrecv_proc];
  status = new MPI_Status[nrecv_proc];

  // work1 = # of datums I send to each proc, including self
  // nsend_proc = # of procs I send messages to, not including self

  for (i = 0; i < nprocs; i++) work1[i] = 0;
  for (i = 0; i < n; i++) work1[proclist[i]]++;

  nsend_proc = 0;
  for (i = 0; i < nprocs; i++)
    if (work1[i]) nsend_proc++;
  if (work1[me]) nsend_proc--;

  // allocate send and self arrays

  proc_send = new int[nsend_proc];
  num_send = new int[nsend_proc];
  index_send = new int[n-work1[me]];
  index_self = new int[work1[me]];

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset work1 to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) {
      num_self = work1[iproc];
      work1[iproc] = 0;
    } else if (work1[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = work1[iproc];
      work1[iproc] = isend;
      isend++;
    }
  }

  // work2 = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // index_self = list of which datums to copy to self

  work2[0] = 0;
  for (i = 1; i < nsend_proc; i++) work2[i] = work2[i-1] + num_send[i-1];

  m = 0;
  for (i = 0; i < n; i++) {
    iproc = proclist[i];
    if (iproc == me) index_self[m++] = i;
    else {
      isend = work1[iproc];
      index_send[work2[isend]++] = i;
    }
  }

  // tell receivers how much data I send
  // sendmax_proc = largest # of datums I send in a single message

  sendmax_proc = 0;
  for (i = 0; i < nsend_proc; i++) {
    MPI_Send(&num_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax_proc = MAX(sendmax_proc,num_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // num_recv = total size of message each proc sends me
  // nrecvdatum = total size of data I recv

  int nrecvdatum = 0;
  for (i = 0; i < nrecv_proc; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvdatum += num_recv[i];
  }
  nrecvdatum += num_self;

  // sort proc_recv and num_recv by proc ID if requested
  // useful for debugging to insure reproducible ordering of received datums

  if (sortflag) {
    int *order = new int[nrecv_proc];
    int *proc_recv_ordered = new int[nrecv_proc];
    int *num_recv_ordered = new int[nrecv_proc];

    for (i = 0; i < nrecv_proc; i++) order[i] = i;
    proc_recv_copy = proc_recv;
    qsort(order,nrecv_proc,sizeof(int),compare_standalone);

    int j;
    for (i = 0; i < nrecv_proc; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      num_recv_ordered[i] = num_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv_proc*sizeof(int));
    memcpy(num_recv,num_recv_ordered,nrecv_proc*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] num_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   comparison function invoked by qsort()
   accesses static class member proc_recv_copy, set before call to qsort()
------------------------------------------------------------------------- */

int compare_standalone(const void *iptr, const void *jptr)
{
  int i = *((int *) iptr);
  int j = *((int *) jptr);
  int *proc_recv = Irregular::proc_recv_copy;
  if (proc_recv[i] < proc_recv[j]) return -1;
  if (proc_recv[i] > proc_recv[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   communicate datums via PlanData
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums (including copied from me)
------------------------------------------------------------------------- */

void Irregular::exchange_data(char *sendbuf, int nbytes, char *recvbuf)
{
  int i,m,n,offset,count;

  // post all receives, starting after self copies

  offset = num_self*nbytes;
  for (int irecv = 0; irecv < nrecv_proc; irecv++) {
    MPI_Irecv(&recvbuf[offset],num_recv[irecv]*nbytes,MPI_CHAR,
              proc_recv[irecv],0,world,&request[irecv]);
    offset += num_recv[irecv]*nbytes;
  }

  // reallocate buf for largest send if necessary

  if (sendmax_proc*nbytes > maxbuf) {
    memory->destroy(buf);
    maxbuf = sendmax_proc*nbytes;
    memory->create(buf,maxbuf,"irregular:buf");
  }

  // send each message
  // pack buf with list of datums
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < nsend_proc; isend++) {
    count = num_send[isend];
    for (i = 0; i < count; i++) {
      m = index_send[n++];
      memcpy(&buf[i*nbytes],&sendbuf[m*nbytes],nbytes);
    }
    MPI_Send(buf,count*nbytes,MPI_CHAR,proc_send[isend],0,world);
  }

  // copy datums to self, put at beginning of recvbuf

  for (i = 0; i < num_self; i++) {
    m = index_self[i];
    memcpy(&recvbuf[i*nbytes],&sendbuf[m*nbytes],nbytes);
  }

  // wait on all incoming messages

  if (nrecv_proc) MPI_Waitall(nrecv_proc,request,status);
}

/* ----------------------------------------------------------------------
   destroy vectors in communication plan for datums
------------------------------------------------------------------------- */

void Irregular::destroy_data()
{
  delete [] proc_send;
  delete [] num_send;
  delete [] index_send;
  delete [] proc_recv;
  delete [] num_recv;
  delete [] index_self;
  delete [] request;
  delete [] status;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void Irregular::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void Irregular::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Irregular::memory_usage()
{
  bigint bytes = 0;
  bytes += maxsend*sizeof(double);   // buf_send
  bytes += maxrecv*sizeof(double);   // buf_recv
  bytes += maxdbuf*sizeof(double);   // dbuf
  bytes += maxbuf;                   // buf
  bytes += 2*maxlocal*sizeof(int);   // mproclist,msizes
  bytes += 2*nprocs*sizeof(int);     // work1,work2
  return bytes;
}
